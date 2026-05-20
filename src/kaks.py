"""Compatibility wrapper for pgkit.commands.kaks."""

from pgkit.commands.kaks import *  # noqa: F401,F403


if __name__ == "__main__":
    from pgkit.commands.kaks import main

    main()
