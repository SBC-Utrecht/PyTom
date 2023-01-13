"""
Export a Conda environment with --from-history, but also append
Pip-installed dependencies
Exports only manually-installed dependencies, excluding build versions, but
including Pip-installed dependencies.
Lots of issues requesting this functionality in the Conda issue tracker, no
sign of progress (as of March 2020).
"""
import subprocess
import sys
import yaml
import argparse


def export_env(history_only=False, include_builds=False):
    """ Capture `conda env export` output """
    cmd = ['conda', 'env', 'export']
    if history_only:
        cmd.append('--from-history')
        if include_builds:
            raise ValueError('Cannot include build versions with "from history" mode')
    if not include_builds:
        cmd.append('--no-builds')
    cp = subprocess.run(cmd, stdout=subprocess.PIPE)
    try:
        cp.check_returncode()
    except Exception as e:
        raise e
    else:
        return yaml.safe_load(cp.stdout)


def _is_history_dep(d, history_deps):
    if not isinstance(d, str):
        return False
    d_prefix = d.split('=')[0]
    return any([d_prefix == dep.split('=')[0] for dep in history_deps if not isinstance(dep, dict)])


def _get_pip_deps(full_deps):
    for dep in full_deps:
        if isinstance(dep, dict) and 'pip' in dep:
            return dep


def _combine_env_data(env_data_full, env_data_hist):
    deps_full = env_data_full['dependencies']
    deps_hist = (env_data_hist['dependencies'])
    deps = [dep for dep in deps_full if _is_history_dep(dep, deps_hist)]
    pip_deps = _get_pip_deps(deps_full)

    env_data = {'name': env_data_full['name']}
    env_data.update(env_data_full)
    env_data['dependencies'] = deps
    env_data['dependencies'].append(pip_deps)

    return env_data


def main(from_history=False, include_builds=False):
    env_data_full = export_env(include_builds=include_builds)
    env_data_hist = export_env(history_only=from_history, include_builds=include_builds)
    env_data = _combine_env_data(env_data_full, env_data_hist)
    yaml.dump(env_data, sys.stdout, sort_keys=False)
    print('Warning: this output might contain packages installed from non-public sources, e.g. a Git repository. '
          'You should review and test the output to make sure it works with `conda env create -f`, '
          'and make changes as required.\n'
          'For example, `conda-env-export` itself is not currently uploaded to PyPI, and it must be removed from '
          'the output file, or else `conda create -f` will fail.', file=sys.stderr)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--from-history', action='store_true',
                        help='Whether to only include packages specified during installation history.')
    parser.add_argument('--include-builds', action='store_true',
                        help='Whether to include builds in addition to version.')
    args = parser.parse_args()
    main(args.from_history, args.include_builds)

