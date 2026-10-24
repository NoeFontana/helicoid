import dataclasses
import json
import tempfile
import unittest
from pathlib import Path

from gen import corpus, manifest, memo
from gen.memo import GEN
from gen.precision import setup
from gen.registry import FUNCTIONS

TASK = ("coeff_k", "theta:exact0")


def fresh_value():
    spec = FUNCTIONS[TASK[0]]
    stratum = next(s for s in spec.strata if s.name == TASK[1])
    return json.loads(json.dumps(list(corpus._build_stratum(spec, stratum))))


class MemoTest(unittest.TestCase):
    def setUp(self):
        setup()
        self.dir = tempfile.TemporaryDirectory()
        self.cache = memo.Memo(self.dir.name)

    def tearDown(self):
        self.dir.cleanup()

    def entry(self) -> Path:
        (path,) = Path(self.dir.name).glob("*.json")
        return path

    def test_a_replay_is_the_fresh_value(self):
        first, how = corpus._registered(*TASK, self.cache)
        self.assertEqual(how, "fresh")
        again, how = corpus._registered(*TASK, self.cache)
        self.assertEqual(how, "hit")
        self.assertEqual(json.loads(json.dumps(first)), again)
        self.assertEqual(again, fresh_value())

    def test_a_task_depends_on_what_it_ran_and_imports_and_not_on_the_registry(self):
        corpus._registered(*TASK, self.cache)
        deps = set(json.loads(self.entry().read_text())["deps"])
        # What it ran, and a module it reads constants of without running it.
        self.assertTrue({"corpus.py", "coeff.py", "strata.py", "__init__.py"} <= deps, deps)
        # The registry names it through the spec, which is in the key; nothing it never imports.
        self.assertNotIn("registry.py", deps)
        self.assertNotIn("sim3.py", deps)

    def test_a_changed_dependency_runs_the_task_again(self):
        corpus._registered(*TASK, self.cache)
        entry = json.loads(self.entry().read_text())
        entry["deps"]["coeff.py"] = "0" * 64
        self.entry().write_text(json.dumps(entry))
        _, how = corpus._registered(*TASK, self.cache)
        self.assertEqual(how, "fresh")

    def test_a_changed_spec_runs_the_task_again(self):
        compute = lambda: [1]
        memo.run(self.cache, ("t", "s"), {"count": 1}, compute)
        self.assertEqual(memo.run(self.cache, ("t", "s"), {"count": 1}, compute)[1], "hit")
        self.assertEqual(memo.run(self.cache, ("t", "s"), {"count": 2}, compute)[1], "fresh")

    def test_an_audit_fails_a_replay_its_fresh_run_does_not_give(self):
        corpus._registered(*TASK, self.cache)
        entry = json.loads(self.entry().read_text())
        entry["value"][1] += 1
        self.entry().write_text(json.dumps(entry))
        audited = dataclasses.replace(self.cache, audit=1.0)
        with self.assertRaises(memo.MemoError):
            corpus._registered(*TASK, audited)

    def test_a_fresh_memo_runs_every_task_and_writes_it(self):
        corpus._registered(*TASK, self.cache)
        entry = json.loads(self.entry().read_text())
        entry["value"][1] += 1
        self.entry().write_text(json.dumps(entry))
        _, how = corpus._registered(*TASK, dataclasses.replace(self.cache, replay=False))
        self.assertEqual(how, "fresh")
        self.assertEqual(json.loads(self.entry().read_text())["value"], fresh_value())

    def test_an_honest_audit_passes(self):
        corpus._registered(*TASK, self.cache)
        _, how = corpus._registered(*TASK, dataclasses.replace(self.cache, audit=1.0))
        self.assertEqual(how, "audited")


class DescribeTest(unittest.TestCase):
    def test_code_closure_and_fields_move_the_description(self):
        def make(c):
            return lambda x: x + c

        self.assertEqual(memo.describe(make(1)), memo.describe(make(1)))
        self.assertNotEqual(memo.describe(make(1)), memo.describe(make(2)))
        self.assertNotEqual(memo.describe(lambda x: x + 1), memo.describe(lambda x: x + 2))
        spec = FUNCTIONS[TASK[0]]
        stratum = spec.strata[0]
        self.assertNotEqual(
            memo.describe(stratum), memo.describe(dataclasses.replace(stratum, count=1))
        )

    def test_every_registered_task_is_described(self):
        for spec in FUNCTIONS.values():
            base = dataclasses.replace(spec, strata=())
            for stratum in spec.strata:
                memo.describe((base, stratum))
        for builder, _ in corpus.SPECIAL.values():
            memo.describe(builder)

    def test_what_it_cannot_describe_is_not_memoised(self):
        with tempfile.TemporaryDirectory() as d:
            opaque = object()
            cache = memo.Memo(d)
            self.assertEqual(memo.run(cache, ("t", "s"), opaque, lambda: 1), (1, "fresh"))
            self.assertEqual(memo.run(cache, ("t", "s"), opaque, lambda: 1), (1, "fresh"))


class DependenciesTest(unittest.TestCase):
    def test_an_import_of_the_registry_is_not_a_dependency(self):
        self.assertNotIn(GEN / "registry.py", memo.dependencies({str(GEN / "corpus.py")}))
        # Its own code running is.
        self.assertIn(GEN / "registry.py", memo.dependencies({str(GEN / "registry.py")}))

    def test_a_module_brings_what_it_imports(self):
        deps = memo.dependencies({str(GEN / "real.py")})
        self.assertIn(GEN / "strata.py", deps)
        self.assertIn(GEN / "precision.py", deps)

    def test_files_outside_gen_are_not_dependencies(self):
        self.assertEqual(memo.dependencies({"/usr/lib/python3/json/__init__.py"}), [GEN / "__init__.py"])



class VerifyTest(unittest.TestCase):
    COMMITTED = GEN.parent.parent / "corpus"

    def copy(self, d: str) -> Path:
        out = Path(d)
        for f in self.COMMITTED.iterdir():
            (out / f.name).write_bytes(f.read_bytes())
        return out

    def test_the_committed_corpus_is_what_its_manifest_says(self):
        self.assertEqual(manifest.verify(self.COMMITTED), [])

    def test_a_changed_missing_or_stray_file_is_named(self):
        with tempfile.TemporaryDirectory() as d:
            out = self.copy(d)
            (out / "so3_log.jsonl").write_bytes(b"{}\n")
            (out / "coeff_k.jsonl").unlink()
            (out / "stray.jsonl").write_text("")
            self.assertEqual(
                manifest.verify(out),
                [
                    "coeff_k.jsonl is listed and missing",
                    "stray.jsonl is present and not listed",
                    "so3_log.jsonl does not have its listed sha256",
                ],
            )

    def test_a_manifest_of_another_generator_is_named(self):
        with tempfile.TemporaryDirectory() as d:
            out = self.copy(d)
            m = json.loads((out / "MANIFEST.json").read_text())
            m["generator"] = "sha256:" + "0" * 64
            (out / "MANIFEST.json").write_text(json.dumps(m))
            (problem,) = manifest.verify(out)
            self.assertTrue(problem.startswith("manifest `generator` is"), problem)

    def test_the_key_moves_with_the_manifest(self):
        with tempfile.TemporaryDirectory() as d:
            out = self.copy(d)
            before = manifest.check_key(out)
            self.assertEqual(before, manifest.check_key(self.COMMITTED))
            (out / "MANIFEST.json").write_bytes((out / "MANIFEST.json").read_bytes() + b" ")
            self.assertNotEqual(manifest.check_key(out), before)


if __name__ == "__main__":
    unittest.main()
