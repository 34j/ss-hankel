import subprocess
import sys

import pytest


@pytest.mark.skipif(
    sys.platform == "win32",
    reason="cyclopts & Windows does not support this test",
)
def test_can_run_as_python_module():
    """Run the CLI as a Python module."""
    result = subprocess.run(  # noqa: S603
        [sys.executable, "-m", "ss_hankel", "--help"],
        check=True,
        capture_output=True,
    )
    assert result.returncode == 0
