## Setting up the environment
Some things need to be prepared if you want to work with the environment of this project.

First open a new terminal.
Then enter the project folder and activate the environment:
```bash
cd cd /data/p_02989/shared_workspace/
${USER}_wd/bat/start.sh results_${USER}   # activates the working environment, make sure your are in that exact directory!
git checkout ${USER}_wd             #check if you are in your local branch

#if you know there were any changes:
git fetch               #check if you missed any changes, alternativly: 
git pull         

cd ${USER}_wd                             # enters the working directory
```
Now open your editor of choice eg Visual studio Code or Pycharm+ or any other editor you like.
```bash
code& #opens Visual Studio Code with the terminal
#alternativly:
pycharm+& #opens pycharm
```
(Note the "&": This guarantees that you can still open any other programm through that same terminal with the same environment.)

## To-Do's when you close the project for today
If you changed anything in the code, don't forget to commit and push your changes. You can do this either in your editor or in your terminal.
For using the terminal use these commands:
```bash
git checkout ${USER}_wd             #check if you are in your local branch
git commit      #commit your changes
git push        #push your changes 
```
For more information, check out [start_with_git.md](/doc/md/start_with_git.md)

[back to README](../../README.md "back to README")