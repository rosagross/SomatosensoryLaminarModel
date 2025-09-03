# How to use git

Git helps us to keep track of all changes and work simultaniously on the same project. the here used Gitlab server is (GWDG): https://gitlab.gwdg.de/r.grossmann/SomatosensoryLaminarModel.git
For initially using Git, follow this recipe in the terminal:
```bash
cd /data/p_02989/Modelling/ # enter the project directory
```

Before you start, make sure you access to gitlab.gwdg.de to this repository: https://gitlab.gwdg.de/r.grossmann/SomatosensoryLaminarModel.git
(An ssh key should be stored in gitLab, otherwise gitlab may ask for passwords whenever so the repository can be cloned. A detailed description can be found here: [Use SSH keys to communicate with GitLab](https://docs.gitlab.com/ee/user/ssh.html "Use SSH keys to communicate with GitLab"))

In a first step, you clone the entire repository into an empty directory to establish your own working directory.
```bash
echo $USER               # checks content of variable USER, this should be username in the institute, eg grossmannr, heurkens, etc 
#export USER=yourname    #  run this to set variable USER if its not your common username
git clone git@gitlab.gwdg.de:r.grossmann/SomatosensoryLaminarModel.git ${USER}_wd   # copying the info from the server into your own working directory
```
This working directory now contains the information of the git server. The current version of the unpacked files is the so-called HEAD of the default branch. Within your working directory, you find subdir .git - please do not delete it, as git needs it for administrative purposes. please work within your own local branch und merge the changes you like to suggest later with the branch main or master. The name of the branch is arbitrary. We suggest to use your username from the institute:
```bash
cd ${USER}_wd 
git branch ${USER}_wd        # introduce a new local branch
git checkout ${USER}_wd      # switch into this branch
```

### Get to your working environment

```bash
cd /data/p_02989/Modelling/
cd ${USER}_wd                             # enters the working directory
```

## Useful codes regarding Gitlab
If you create new files, you need to stage them:
```bash
git add <your-new-file>
```
At certain time points you should commit the changes within your working directory 
to the local git management. 
```bash
git commit -a -m 'some reasonable comment' # -a commits all staged changes
```
To inform the server about the local commits, you need to push the information there:
```bash
git push
```

Check for new information at the server:
```bash
git fetch --all                   # retrieves information about all registered remote directories and their branches
```
Checkin you current changes locally:
```bash
git stash                         # preserve local changes in the stack
```
or remove the local changes
```bash
git checkout .     # or
git restore .
```
Integrate the current main repo into your branch (works from the older branch to the newer). 
Make sure you committed/pushed all your changes beforehand.
```bash
git merge main # provided main is newer than your current branch
git push # inform the server about the merge
```
Otherwise checkout main and merge it into your branch:
```bash
git checkout main
git pull
git merge <your-branch> # provided your current branch is newer than main
git push # inform the server about the merge
git checkout <your-branch> # switch back to your branch 
```

[back to README](../../README.md "back to README")
