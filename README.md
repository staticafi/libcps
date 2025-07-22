# **CPS**

**CPS** is an abreviation for **C**overage **P**problem **S**olver.

## **Dependencies**

**libutility**: https://github.com/staticafi/libutility

**Eigen3**: https://gitlab.com/libeigen/eigen/tree/master
- We recommend to install it via `vcpkg`.
    - If you use VS Code, then:
        - Copy the file `setup/VSCode/vcpkg.json` to root
        of your entire project.
        - Copy the file `setup/VSCode/setting.json` to .vscode folder under the root folder of your entire project.
        Also, replace `<vcpkg-install-dir>` by valid installation directory of `vcpkg`.

- In VS Code you can setup `Eigen`'s pretty printer for `GDB`:
    - Add this under the `"setupCommands"` section in your
    `launch.json` file (with `<your-name>` replaced, of course):
        ```
        {
            "description": "Source .gdbinit manually for Eigen pretty-printers",
            "text": "source /home/<your-name>/.gdbinit"
        }
        ```
    - Open (or create, if not present) `.gdbinit` file in your home directory
    and include there this text:
        ```
        python
        import sys
        sys.path.insert(0, "/path/to/eigen/printer/")
        from printers import register_eigen_printers
        register_eigen_printers(None)
        end
        ```
    The path must be the directory into which you downloaded files from here: 
    https://gitlab.com/libeigen/eigen/-/tree/master/debug/gdb