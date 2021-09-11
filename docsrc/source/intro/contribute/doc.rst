How to extend the documentation
================================



The documentation is build using sphinx and deployed to github using github Actions.
To extend the documentation

1. Make changes in the docsrc/source directory

2. Locally build the documentation

    .. code-block::

        cd docsrc
        make html

3. Inspect resulting docsrc/build/html/index.html file

4. When satisfied push changes in docsrc/source to master.

    .. warning::

        Only add the docsrc/source folder to github, do not commit the build folder.

5. A github action will be trigger regularly to deploy your changes to github pages



Notes
-----
- If you have trouble building the documentation, make sure you have all dependencies installed. The fenics_explicit.yml file includes all dependencies for building this documentation with sphinx.

- The documentation is deployed to the gp-pages branch using https://github.com/peaceiris/actions-gh-pages

- Sphinx allows building to many different formats. Try

    .. code-block::

        make latex

External resources
---------------------

- Sphinx documentation: https://www.sphinx-doc.org/en/master/
- Introduction to reStructuredText: https://www.sphinx-doc.org/en/master/usage/restructuredtext/basics.html