from ss_hankel.cli import app


def test_help():
    """The help message includes the CLI name."""
    app(["--help"])


def test_cli():
    app(
        [
            "{{3+Exp[x/2],2+2x+Exp[x/2]},{3+Exp[x/2],-1+x+Exp[x/2]}}",
            "--circle-radius",
            "4",
        ]
    )
