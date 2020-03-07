import nox

@nox.session
def lint(session):
    session.install("flake8")
    session.run("flake8", "./analyzers")
    session.run("flake8", "./experiments")
    session.run("flake8", "./model")
