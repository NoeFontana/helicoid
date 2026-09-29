import hashlib
import json
import tempfile
import unittest
from pathlib import Path
from unittest import mock

from mpmath import mp, mpf

from gen import corpus, manifest, precision
from gen.registry import FUNCTIONS, FunctionSpec
from gen.strata import Stratum

precision.setup()
COMMITTED = manifest.ROOT.parent / "corpus"


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

    def test_ids_are_sequential_and_strata_in_canonical_order(self):
        data, n, _ = corpus.build(FUNCTIONS["coeff_k"], only={"theta:1e-3", "theta:exact0"})
        records = [json.loads(line) for line in data.decode().splitlines()]
        self.assertEqual([r["id"] for r in records], list(range(n)))
        self.assertEqual([r["stratum"] for r in records], ["theta:1e-3"] * 64 + ["theta:exact0"])


class GenerationTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.dirs = []
        for _ in range(2):
            tmp = tempfile.TemporaryDirectory()
            cls.addClassCleanup(tmp.cleanup)
            corpus.write(Path(tmp.name))
            cls.dirs.append(Path(tmp.name))

    def test_generation_is_deterministic(self):
        a, b = self.dirs
        names = sorted(p.name for p in a.iterdir())
        self.assertEqual(names, sorted(p.name for p in b.iterdir()))
        for name in names:
            self.assertEqual((a / name).read_bytes(), (b / name).read_bytes(), name)

    def test_regeneration_reproduces_the_committed_corpus(self):
        """Behaviour, not provenance: every byte but the manifest's `generator` field."""
        regenerated = self.dirs[0]
        self.assertEqual(
            sorted(p.name for p in regenerated.iterdir()),
            sorted(p.name for p in COMMITTED.iterdir()),
        )
        for path in regenerated.glob("*.jsonl"):
            self.assertEqual(path.read_bytes(), (COMMITTED / path.name).read_bytes(), path.name)
        fresh, committed = (
            json.loads((d / "MANIFEST.json").read_text()) for d in (regenerated, COMMITTED)
        )
        fresh.pop("generator"), committed.pop("generator")
        self.assertEqual(fresh, committed)

    def test_write_removes_files_no_function_id_owns(self):
        with tempfile.TemporaryDirectory() as d:
            out = Path(d)
            (out / "renamed_away.jsonl").write_text("stale\n")
            (out / "notes.txt").write_text("not ours\n")
            names = set(corpus.write(out))
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
        m = json.loads(manifest.render({"f.jsonl": (b"x\n", 1, 1)}))
        self.assertEqual(
            m["files"]["f.jsonl"],
            {"records": 1, "rechecked": 1, "sha256": hashlib.sha256(b"x\n").hexdigest()},
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

    def test_committed_manifest_names_this_generator(self):
        """Provenance only: any edit to the sources fails it; it says nothing about behaviour."""
        m = json.loads((COMMITTED / "MANIFEST.json").read_text())
        self.assertEqual(
            m["generator"], manifest.generator_identity(), "generator changed: run `just corpus`"
        )


if __name__ == "__main__":
    unittest.main()
