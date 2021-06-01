from invoke import task, call

@task(aliases=['del'])

def delete(c, extra=''):
    patterns = ['*.so']
    if extra:
        patterns.append(extra)
    for pattern in patterns:
        c.run("rm -rf {}".format(pattern))
    print("Done deleting")

@task(delete)
def build(c):
    c.run("python3.8.5 setup.py build_ext --inplace")
    print("Done building")