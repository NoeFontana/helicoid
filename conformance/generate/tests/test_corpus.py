import dataclasses
import hashlib
import json
import tempfile
import unittest
from pathlib import Path
from unittest import mock

from mpmath import mp, mpf

from gen import corpus, fmt, manifest, precision
from gen.manifest import Built
from gen.registry import FUNCTIONS, FunctionSpec
from gen.strata import SCALAR_THETA_STRATA, Stratum

precision.setup()
COMMITTED = manifest.ROOT.parent / "corpus"
ASKED_FOR_F32 = ("solve_cubic", "eig3", "chol_", "quat_", "real_", "mat2_")  # 0056, 0061
SUBSET_STRATUM = next(s for s in SCALAR_THETA_STRATA if s.name == "theta:exact0")


def noisy_spec(noise: str, n: int = 1, bad=range(1)) -> FunctionSpec:
    """`v` is 1/3, except that at 150 digits inputs in `bad` are off by `noise`."""

    def evaluate(inp: dict) -> dict:
        off = mpf(noise) if mp.dps > 130 and inp["i"] in bad else 0
        return {"v": mpf(1) / 3 + off}

    stratum = Stratum("s", n, lambda rng: [float(i) for i in range(n)])
    return FunctionSpec("t", (stratum,), lambda st: [{"i": i} for i in range(n)], evaluate)


class RecheckTest(unittest.TestCase):
    def test_disagreement_fails_generation(self):
        ok, bad = noisy_spec("1e-45"), noisy_spec("1e-30")
        corpus.recheck(ok, {"i": 0}, {"v": mpf(1) / 3})
        with self.assertRaises(corpus.RecheckError):
            corpus.recheck(bad, {"i": 0}, {"v": mpf(1) / 3})

    def test_shape_disagreement_fails_generation(self):
        spec = FunctionSpec("t", (), lambda s: [], lambda i: {"v": [mpf(1), mpf(2)]})
        with self.assertRaises(corpus.RecheckError):
            corpus.recheck(spec, {}, {"v": [mpf(1)]})
        with self.assertRaises(corpus.RecheckError):
            corpus.recheck(spec, {}, {"w": [mpf(1), mpf(2)]})

    def test_indices_cover_every_stratum_and_one_percent(self):
        for n in (1, 2, 64, 100, 101, 801):
            picks = corpus.recheck_indices(n)
            self.assertTrue(picks and all(0 <= i < n for i in picks), n)
            self.assertGreaterEqual(len(picks), n / 100)


class BuildTest(unittest.TestCase):
    def test_golden_record(self):
        data, n, rechecked = corpus.build(FUNCTIONS["coeff_k"], only={"theta:exact0"})
        self.assertEqual((n, rechecked), (1, 1))
        self.assertEqual(
            data,
            b'{"id":0,"in":{"theta":"0x0.0p+0"},"out":{"d_branch":"-2.08333333333333333333333333333e-2",'
            b'"value":"5.00000000000000000000000000000e-1"},"stratum":"theta:exact0"}\n',
        )

    def test_golden_records_of_the_f32_strata(self):
        """Each `@f32` record is the function at a binary32 input; `q:w0@f32` rounds n = 1e-3."""
        data, n, _ = corpus.build(FUNCTIONS["coeff_cos_half"], only={"theta:exact0@f32"})
        self.assertEqual(
            data,
            b'{"id":0,"in":{"theta":"0x0.0p+0"},"out":{"d_branch":"-1.25000000000000000000000000000e-1",'
            b'"value":"1.00000000000000000000000000000e0"},"stratum":"theta:exact0@f32"}\n',
        )
        data, n, _ = corpus.build(FUNCTIONS["coeff_r"], only={"q:w0@f32"})
        second = json.loads(data.decode().splitlines()[1])
        self.assertEqual((n, second["in"]["n"]), (3, "0x1.0624de0000000p-10"))  # f32(1e-3)
        with mp.workdps(60):  # r = pi / n at the rounded n, not at 1e-3
            n32, r = mpf(float.fromhex(second["in"]["n"])), mpf(second["out"]["value"])
            self.assertLess(abs(r - mp.pi / n32), mpf(10) ** -28 * mp.pi / n32)
            self.assertGreater(abs(r - mp.pi / mpf("1e-3")), mpf(10) ** -4)

    def test_an_f32_stratum_with_a_non_binary32_input_fails_generation(self):
        stratum = Stratum("s@f32", 1, lambda rng: [0.1], f32=True)
        one = lambda inp: {"v": mpf(1)}
        for x in (0.1, [0.5, 1e-310]):
            spec = FunctionSpec("t", (stratum,), lambda st, x=x: [{"x": x}], one)
            with self.assertRaises(corpus.BinaryError, msg=x):
                corpus.build(spec)
        exact = FunctionSpec("t", (stratum,), lambda st: [{"x": [0.5, 2.0**-149]}], one)
        self.assertEqual(corpus.build(exact)[1], 1)

    def test_golden_records_of_r_at_w_zero(self):
        data, n, _ = corpus.build(FUNCTIONS["coeff_r"], only={"q:w0"})
        self.assertEqual(n, 3)
        first, *rest = data.splitlines(keepends=True)
        self.assertEqual(
            first,
            b'{"id":0,"in":{"n":"0x1.0000000000000p+0","w":"0x0.0p+0"},"out":{"d_branch":'
            b'"-1.57079632679489661923132169164e0","value":"3.14159265358979323846264338328e0"},'
            b'"stratum":"q:w0"}\n',
        )
        for line, norm in zip(rest, (1e-3, 1e3), strict=True):  # r = pi / n, r' = -pi / (2 n^3)
            record = json.loads(line)
            self.assertEqual(float.fromhex(record["in"]["n"]), norm)
            with mp.workdps(60):
                want = (mp.pi / mpf(norm), -mp.pi / (2 * mpf(norm) ** 3))
                got = (mpf(record["out"]["value"]), mpf(record["out"]["d_branch"]))
            for g, w in zip(got, want, strict=True):
                self.assertLess(abs(g - w), mpf(10) ** -29 * abs(w), (norm, g, w))

    def test_golden_records_of_the_matrix_encoding(self):
        """A `Mat` is a column-major array with a sibling `shape`, as an output and as an input."""
        one, zero = "1.00000000000000000000000000000e0", "0.00000000000000000000000000000e0"
        eye = ",".join(f'"{one if i % 4 == 0 else zero}"' for i in range(9))
        data, _, _ = corpus.build(FUNCTIONS["so3_jr"], only={"theta:exact0"})
        self.assertEqual(
            data,
            (
                '{"id":0,"in":{"phi":["0x0.0p+0","0x0.0p+0","0x0.0p+0"]},"out":{"J":[%s],'
                '"shape":[3,3]},"stratum":"theta:exact0"}\n' % eye
            ).encode(),
        )
        data, _, _ = corpus.build(FUNCTIONS["so3_from_matrix"], only={"theta:exact0"})
        (first,) = data.splitlines(keepends=True)
        hexes = ",".join(f'"{float(i % 4 == 0).hex()}"' for i in range(9))
        self.assertEqual(
            first,
            (
                '{"id":0,"in":{"R":[%s],"shape":[3,3]},"out":{"q":["%s","%s","%s","%s"]},'
                '"stratum":"theta:exact0"}\n' % (hexes, one, zero, zero, zero)
            ).encode(),
        )

    def test_two_matrices_in_one_record_part_are_refused(self):
        mat = fmt.Mat((1, 1), [1.0])
        with self.assertRaises(ValueError):
            corpus._encode({"a": mat, "b": mat}, fmt.hex_float)

    def test_build_runs_the_spec_check_on_every_record(self):
        seen = []
        spec = dataclasses.replace(noisy_spec("0", 3), check=lambda inp, out: seen.append(inp["i"]))
        corpus.build(spec)
        self.assertEqual(seen, [0, 1, 2])

        def refuse(inp, out):
            raise RuntimeError("disagreement")

        with self.assertRaises(RuntimeError):
            corpus.build(dataclasses.replace(noisy_spec("0", 3), check=refuse))

    def test_build_rechecks_and_counts_what_it_checked(self):
        # 250 records: three rechecked, at the indices of `recheck_indices`.
        picks = sorted(corpus.recheck_indices(250))
        self.assertEqual(len(picks), 3)
        _, n, rechecked = corpus.build(noisy_spec("1e-45", 250, bad=picks))
        self.assertEqual((n, rechecked), (250, 3))
        with self.assertRaises(corpus.RecheckError):
            corpus.build(noisy_spec("1e-35", 250, bad=picks[-1:]))
        with self.assertRaises(corpus.RecheckError):
            corpus.build(noisy_spec("1e-35"))

    def test_the_f32_strata_come_after_every_binary64_one_in_the_same_order(self):
        for name in ("coeff_k", "coeff_a", "coeff_b", "coeff_c", "coeff_d", "coeff_e"):
            strata = [s.name for s in FUNCTIONS[name].strata]
            half = len(strata) // 2
            self.assertEqual(strata[half:], [f"{s}@f32" for s in strata[:half]], name)
            self.assertEqual(strata[:half], [s.name for s in SCALAR_THETA_STRATA], name)
        strata = [s.name for s in FUNCTIONS["coeff_r"].strata]
        self.assertEqual(strata[:29], [*(s.name for s in SCALAR_THETA_STRATA), "q:w0"])
        self.assertEqual(strata[29:], [*(f"{s.name}@f32" for s in SCALAR_THETA_STRATA), "q:w0@f32"])
        other = [n for n in FUNCTIONS if not n.startswith(("coeff_", *ASKED_FOR_F32))]
        for name in other:  # no other id has an @f32 stratum until a record asks for one
            self.assertFalse(any(s.f32 for s in FUNCTIONS[name].strata), name)

    def test_the_0056_ids_have_their_f32_strata_after_every_binary64_one(self):
        """An `@f32` stratum is named for its binary64 namesake where one exists (0056 decision 2)
        and follows the binary64 order; the others (`x:1e37@f32`) are the binary32 range's own."""
        names = [n for n in FUNCTIONS if n.startswith(ASKED_FOR_F32)]
        self.assertEqual(len(names), 13 + 1)  # and 0061's `mat2_inverse_adj`
        for name in names:
            strata = FUNCTIONS[name].strata
            flags = [s.f32 for s in strata]
            self.assertEqual(flags, sorted(flags), name)  # every binary64 stratum first
            half = [s.name for s in strata if not s.f32]
            twins = [s.name.removesuffix("@f32") for s in strata if s.f32]
            self.assertTrue(all(s.name.endswith("@f32") == s.f32 for s in strata), name)
            shared = [t for t in twins if t in half]
            self.assertEqual(shared, sorted(shared, key=half.index), name)
            for s in strata:  # a twin reads its namesake's stream exactly when it has one
                if s.f32 and s.stream_of is not None:
                    self.assertEqual(s.stream_of, s.name.removesuffix("@f32"), (name, s.name))

    def test_ids_are_sequential_and_strata_in_canonical_order(self):
        data, n, _ = corpus.build(FUNCTIONS["coeff_k"], only={"theta:1e-3", "theta:exact0"})
        records = [json.loads(line) for line in data.decode().splitlines()]
        self.assertEqual([r["id"] for r in records], list(range(n)))
        self.assertEqual([r["stratum"] for r in records], ["theta:1e-3"] * 64 + ["theta:exact0"])


# The full corpus is regenerated and compared by `just corpus-check`; the tests regenerate these.
SUBSET = {"theta:exact0", "theta:1e-6", "q:w0", "theta:exact0@f32", "theta:1e-6@f32", "q:w0@f32"}
# and one cheap stratum of each 0056 family, a shared-stream and an own-stream @f32 among them
SUBSET |= {"cubic:one-real", "eig:triple@f32", "renorm:eta-edge@f32", "nd:generic@f32", "nd:wide"}
SUBSET |= {"x:1e-37@f32", "yx:ratio-1e+8", "chol:cond-1e-2@f32"}


def subset_registry() -> dict:
    keep = lambda spec: tuple(s for s in spec.strata if s.name in SUBSET)
    return {name: dataclasses.replace(spec, strata=keep(spec)) for name, spec in FUNCTIONS.items()}


def records(data: bytes) -> list[dict]:
    """The records of a file without their ids, which count from the start of the file."""
    return [{k: v for k, v in json.loads(line).items() if k != "id"} for line in data.splitlines()]


class GenerationTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        with mock.patch.object(corpus, "FUNCTIONS", subset_registry()):
            cls.serial, cls.pooled = corpus.build_all(jobs=1), corpus.build_all(jobs=2)

    def test_generation_is_deterministic_whatever_the_worker_count(self):
        self.assertEqual(self.serial, self.pooled)

    def test_regeneration_reproduces_the_committed_records(self):
        self.assertEqual(sorted(self.serial), sorted(FUNCTIONS))
        for name, (data, n, _) in self.serial.items():
            committed = records((COMMITTED / f"{name}.jsonl").read_bytes())
            want = [r for r in committed if r["stratum"] in SUBSET]
            self.assertEqual(len(want), n, name)
            self.assertEqual(records(data), want, name)

    def test_write_removes_files_no_function_id_owns(self):
        tiny = {"coeff_k": dataclasses.replace(FUNCTIONS["coeff_k"], strata=(SUBSET_STRATUM,))}
        with tempfile.TemporaryDirectory() as d, mock.patch.object(corpus, "FUNCTIONS", tiny):
            out = Path(d)
            (out / "renamed_away.jsonl").write_text("stale\n")
            (out / "notes.txt").write_text("not ours\n")
            names = set(corpus.write(out, jobs=1))
            self.assertEqual(
                sorted(p.name for p in out.iterdir()),
                sorted([*names, "MANIFEST.json", "notes.txt"]),
            )


class ManifestTest(unittest.TestCase):
    def test_identity_covers_sources_and_nothing_else(self):
        with tempfile.TemporaryDirectory() as d:
            root = Path(d)
            (root / "gen").mkdir()
            for name in ("pyproject.toml", "uv.lock", ".python-version", "gen/a.py"):
                (root / name).write_text(name)
            base = manifest.generator_identity(root)
            (root / "gen/notes.txt").write_text("ignored")
            self.assertEqual(manifest.generator_identity(root), base)
            for name in ("uv.lock", "gen/a.py", ".python-version"):
                (root / name).write_text("changed")
                self.assertNotEqual(manifest.generator_identity(root), base, name)
                (root / name).write_text(name)
            (root / "gen/b.py").write_text("")
            self.assertNotEqual(manifest.generator_identity(root), base)

    def test_identity_is_sorted_framed_and_independent_of_directory_order(self):
        names = (".python-version", "gen/a.py", "gen/sub/c.py", "pyproject.toml", "uv.lock")
        with tempfile.TemporaryDirectory() as d:
            root = Path(d)
            (root / "gen/sub").mkdir(parents=True)
            for name in reversed(names):
                (root / name).write_text(name)
            expected = hashlib.sha256()
            for name in names:
                data = (root / name).read_bytes()
                expected.update(f"{name}\0{len(data)}\0".encode() + data + b"\n")
            self.assertEqual(manifest.generator_identity(root), "sha256:" + expected.hexdigest())
            listed = Path.rglob
            with mock.patch.object(Path, "rglob", lambda p, pat: reversed(sorted(listed(p, pat)))):
                self.assertEqual(
                    manifest.generator_identity(root), "sha256:" + expected.hexdigest()
                )

    def test_render_records_hashes_counts_and_versions(self):
        m = json.loads(
            manifest.render(
                {"f.jsonl": Built(b"x\n", 1, 1), "s.jsonl": Built(b"y\n", 1, 1, kind="series")}
            )
        )
        sha = lambda data: hashlib.sha256(data).hexdigest()
        self.assertEqual(
            m["files"],
            {
                "f.jsonl": {"kind": "corpus", "records": 1, "rechecked": 1, "sha256": sha(b"x\n")},
                "s.jsonl": {"kind": "series", "records": 1, "verified": 1, "sha256": sha(b"y\n")},
            },
        )
        self.assertEqual(
            (m["dps"], m["seed"], m["recheck"]),
            (120, "0x68656c69636f6964", {"digits": 40, "dps": 150}),
        )
        self.assertRegex(m["python"], r"^3\.12$")

    def test_committed_files_match_their_manifest_entries(self):
        m = json.loads((COMMITTED / "MANIFEST.json").read_text())
        self.assertEqual(
            sorted(p.name for p in COMMITTED.iterdir()), sorted([*m["files"], "MANIFEST.json"])
        )
        for name, entry in m["files"].items():
            data = (COMMITTED / name).read_bytes()
            self.assertEqual(hashlib.sha256(data).hexdigest(), entry["sha256"], name)
            self.assertEqual(data.count(b"\n"), entry["records"], name)

    def test_committed_kinds_say_what_the_checked_count_means(self):
        """`rechecked` is records recomputed at 150 digits; `coeff_series` has `verified` instead.

        Three kinds, and only `corpus` holds (in, out, stratum) records a subject is scored on:
        `series` supplies the generated files' literals and `switch-ref` the right-hand side of
        docs/maths/coefficients.md CO.12 (docs/decisions/0039 plan step 0).
        """
        kinds = {"coeff_series.jsonl": "series", "coeff_switch_ref.jsonl": "switch-ref"}
        files = json.loads((COMMITTED / "MANIFEST.json").read_text())["files"]
        for name, entry in files.items():
            kind = kinds.get(name, "corpus")
            self.assertEqual(entry["kind"], kind, name)
            self.assertEqual("verified" in entry, kind == "series", name)
            self.assertEqual("rechecked" in entry, kind != "series", name)
            for line in (COMMITTED / name).read_text().splitlines():
                record = json.loads(line)
                has = {"in", "out", "stratum"} <= record.keys()
                self.assertEqual(has, kind == "corpus", name)

    def test_committed_manifest_names_this_generator(self):
        """Provenance only: any edit to the sources fails it; it says nothing about behaviour."""
        m = json.loads((COMMITTED / "MANIFEST.json").read_text())
        self.assertEqual(
            m["generator"], manifest.generator_identity(), "generator changed: run `just corpus`"
        )


if __name__ == "__main__":
    unittest.main()
