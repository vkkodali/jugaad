# AI Instructions

- Use `uv` for dependency, environment, build, test, and release commands.
- Keep CLI behavior stable unless the user asks for a breaking change.
- Run `uv run pytest`, `uv run ruff check .`, and `uv run ruff format --check .` before finishing code changes.
- Use Conventional Commits for commit messages; semantic-release depends on them.
