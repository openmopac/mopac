# Viewing MOPAC HTML output with JSmol

When a MOPAC input contains the `HTML` keyword, MOPAC writes a web page
(`name.html`) that displays the structure in JSmol. Since MOPAC 23.2.1 that page
loads JSmol from `/jsmol/JSmol.min.js`. Modern browsers do not let a page opened
directly from disk (a `file://` address) load JSmol, so the page has to be
served by a small local web server whose root folder contains a `jsmol` folder.

The scripts here do that for you: open an HTML file with them and the molecule
appears in your default browser. Nothing is copied into your working folders,
no administrator rights are needed, and the server only listens on this
computer (`localhost`).

The MOPAC installer puts these scripts in a `jsmol-launcher` folder in the MOPAC
installation. If you installed MOPAC another way (for example with conda), copy
this folder from <https://github.com/openmopac/mopac/tree/main/jsmol-launcher>. The scripts do not depend on the
MOPAC version.

## One-time setup: put one copy of JSmol above your working folders

Download the latest `Jmol-xx.x.x-binary.zip` from
<https://sourceforge.net/projects/jmol/files/latest/download>, open it, open the
`jsmol.zip` inside it, and copy the `jsmol` folder out into a folder that sits
**above** the folders where you run MOPAC. For example:

```
C:\Users\you\MOPAC\              (or ~/MOPAC on macOS/Linux)
    jsmol\                       <- contains JSmol.min.js, j2s\, ...
    Project1\
        ethanol.html
    Project2\Runs\today\
        test.html
```

Working folders can be at any depth below the folder holding `jsmol`.

If you would rather keep JSmol somewhere else, set the environment variable
`MOPAC_JSMOL_ROOT` to the folder that *contains* `jsmol`. It is used when no
`jsmol` folder is found above the HTML file. Note that the HTML file must then
also be somewhere below `MOPAC_JSMOL_ROOT`, because the server can only serve
files under its root.

## Windows

The scripts are in the `jsmol-launcher\windows` folder (with the MOPAC
installer, `C:\Program Files\MOPAC\jsmol-launcher\windows`).

* Drag an HTML file onto `Open_MOPAC_HTML.cmd`, or
* from a command prompt: `Open_MOPAC_HTML.cmd "C:\Users\you\MOPAC\Project1\ethanol.html"`, or
* right-click a `.html` file, choose "Open with", "Choose another app", browse
  to `Open_MOPAC_HTML.cmd`, and tick "Always use this app".

The first time, a minimized window titled "MOPAC JSmol server" starts; it is
reused for later files under the same root, so they open immediately. Stop it
with `Stop_MOPAC_server.cmd` or by closing that window. It uses only PowerShell,
which is part of Windows; Python is not needed.

## macOS and Linux

```
jsmol-launcher/unix/open-mopac-html.sh ~/MOPAC/Project1/ethanol.html
jsmol-launcher/unix/open-mopac-html.sh --stop
```

Requires `python3`. Options: `--port N` to start at another port (default
8000; the next free port is used if it is taken) and `--no-open` to print the
address without opening a browser.

## How it works

1. Starting from the folder holding the HTML file, the launcher walks up the
   folder tree until it finds a folder containing `jsmol/JSmol.min.js`.
2. It starts a local web server rooted at that folder, unless one for the same
   folder is already running.
3. It opens `http://localhost:<port>/<path of the file below that folder>` in
   the default browser, so `/jsmol/JSmol.min.js` in the page resolves to the one
   shared copy of JSmol.

The server keeps running so that JSmol can load the extra files it fetches on
demand (for example when a menu is first opened).

You can do the same by hand with any static web server, e.g.
`python -m http.server` started in the folder that holds `jsmol`; see
<https://openmopac.net/guides/jsmol/>.

## Older MOPAC output

HTML written by MOPAC 23.2 and earlier loads `../jsmol/JSmol.min.js`, i.e. a
`jsmol` folder next to the folder holding the HTML file. The launcher serves
those pages too, as long as the HTML file is exactly one folder below the
folder holding `jsmol`.
