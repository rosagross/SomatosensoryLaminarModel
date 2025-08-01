# Directory Management
Please note that when working with the project code, you should use aliases for your used directories. 
Here you can look up which directory should contain what kind of files and how you can use them.

## Dir in the code
- EXPDIR 	-->	"/shared_workspace/" 
- RAWDIR 	-->	"/shared_workspace/rawdir" 
- DOCDIR 	-->	"/shared_workspace/doc" 
- RESDIR 	-->	"/shared_workspace/results_$USER/" 
- WDDIR	-->	"/shared_workspace/$USER_wd/" 
- DATADIR	-->	"/shared_workspace/datadir/" 
- PERSDATA --> "/shared_workspace/$USER_wd/data" 

## folder in project
- /shared_workspace/rawdir/
    --> for raw, unedited data
- /shared_workspace/$USER_wd/
     --> for everything thats going to be worked with: code, data, figures etc.
- /shared_workspace/datadir/
    --> for edited but reused data, shared with everyone else
- /shared_workspace/$USER_wd/code
    --> any form of code-files
- /shared_workspace/$USER_wd/code/lib
    --> helper files 
- /shared_workspace/$USER_wd/doc 
    -->any form of document like papers, git-documentation or similar
- /shared_workspace/$USER_wd/data
    --> data you use for yourself, not shared with others
- /shared_workspace/results_$USER/ 
    --> for everything thats not going to be worked with anymore: figures, tables, results etc. a new subfolder for every subject

[back to README](../../README.md "back to README")