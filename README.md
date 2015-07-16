ngs-utils
=========

Collection of notes and scripts related to NGS.



# Make git work on our cluster:
git clone https://github.com/YOUR_USER_NAME/ngs-utils.git

# Make git work via ssh
NB This does not work yet. Probably we have to upgrade our git module?

Look up your git email (go to github.com, click your photo / 
settings / email).

```bash
cd
module load git/1.9.3
ssh-keygen -t rsa -b 4096 -C "your_email@example.com"
mkdir git

git clone git@github.com:YOUR_USER_NAME/ngs-utils.git


cd ngs-utils
git config --global user.name "YOUR NAME"
git config --global user.email "YOUR EMAIL ADDRESS"

eval "$(ssh-agent -s)"     # start ssh-agent in bg
ssh-add ~/.ssh/id_rsa      # add your key to the ssh-agent
```
Now copy the public key (~/.ssh/id_rsa.pub) exactly without adding newlines or whitespace to your 
clipboard!

On github.com, click you photo, goto settings / Add SSH-key, and add the key and test your connection:
```
ssh -T git@github.com
```
