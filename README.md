# Cime

This is the repository for the public cime library (as discussed in the paper). It builds upon the PSE library found under the following repository https://github.com/jku-vds-lab/projection-space-explorer.

## Installation

First clone this repository, preferably using ssh.

```
git clone git@github.com:jku-vds-lab/cime.git
```

After this step you can install the dependencies using

```
npm install
```

Note that PSE has several peer dependencies which also need to be installed (material ui, react, reactdom etc). Either look them up in the package.json and install them, or use a library which automatically manages that.

## Linking PSE

If you want to make changes to PSE and view the changes without having to push to the repo and reinstalling dependencies, the recommended way is to use the npm-link feature. For this to work you first need to clone the PSE repository using

```
git clone git@github.com:jku-vds-lab/projection-space-explorer.git
```

Then navigate to the root folder and create a symlink using

```
npm link
```

After this you can navigate to your cime folder and link to your local PSE version using

```
npm link projection-space-explorer
```

(Note that you need to to this step AFTER the npm install step)

After this step you can run your changes in PSE, then build them using the commands provided in the PSE repository.