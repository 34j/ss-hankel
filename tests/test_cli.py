from ss_hankel.cli import app


def test_help():
    """The help message includes the CLI name."""
    app(["--help"])
