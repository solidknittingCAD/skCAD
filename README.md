# skCAD
![skCAD pipeline showing blocks converted into a knitting pattern, and the resulting physical Stanford bunny](https://textiles.cs.cmu.edu/publications/2026-skcad/teaser.jpg)

A design tool for solid knitting with automatic pattern generation.

* Project page: https://textiles.cs.cmu.edu/publications/2026-skcad/
* Paper: https://doi.org/10.1145/3772318.3791148

## Video Overview
https://youtu.be/4rld9zOaUzU

## Running the code
Go to https://solidknittingCAD.github.io/skCAD/

## Editing Patterns
You can create 3D forms by stacking unlabeled blocks using the commands below.

* `A`: add a block at your current cursor location
* `X`: delete the block under your cursor

After finishing the design, click `Edit/Pattern` to convert the blocks into a solid-knitting pattern.

In Pattern mode, you can highlight a specific loop by clicking `Indicate Loop`.  
Press `<` or `>` on the keyboard to move to the previous or next loop.

Click the settings gear in the top-left corner, and use the checkboxes to show or hide specific layers.

Click `Load` or `Save` to load or save a block assembly as a .body file.  
For example, try loading a Stanford bunny .body file from [here](https://solidknittingCAD.github.io/skCAD/patterns/stanford-bunny.body).

Click `Save Yarn Path` to export the pattern as an .obj file.

## Notes
If you press `<` or `>` immediately after clicking the checkbox, the checkbox will remain selected and the keys won’t work.  
Click anywhere near the top of the page to deselect it.

## Acknowledgment
This code is based on https://github.com/textiles-lab/solid-knitting-ui