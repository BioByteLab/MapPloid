"""Supply controlled GenMap command outputs for region/filter unit tests."""
from contextlib import contextmanager
from pathlib import Path
from unittest.mock import patch
import shutil


@contextmanager
def mock_genmap(scores):
    original_which = shutil.which

    def which(name, *args, **kwargs):
        return "/fixture/genmap" if name == "genmap" else original_which(name, *args, **kwargs)

    def run(command, *args, **kwargs):
        if command[1] == "index":
            index = Path(command[command.index("-I") + 1])
            index.mkdir(parents=True, exist_ok=True)
            (index / "index.sa.val").write_text("fixture")
        elif command[1] == "map":
            target = Path(command[command.index("-O") + 1] + ".bedgraph")
            shutil.copyfile(scores, target)
        else:
            raise AssertionError(command)

    with patch("mapploid.reference_mask.shutil.which", side_effect=which), \
         patch("mapploid.reference_mask.run_command", side_effect=run):
        yield
