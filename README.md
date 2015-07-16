ngs-utils
=========

Collection of notes and scripts related to NGS.



Steps to make git work on our cluster:

# Look up your git email (go to github.com, click your photo / 
settings / email)


cd
module load git/1.9.3
ssh-keygen -t rsa -b 4096 -C "your_email@example.com"
mkdir git

git clone git@github.com:YOUR_USER_NAME/ngs-utils.git
# Or, if you want via https: git clone https://github.com/YOUR_USER_NAME/ngs-utils.git 

cd ngs-utils
git config --global user.name "YOUR NAME"
git config --global user.email "YOUR EMAIL ADDRESS"

eval "$(ssh-agent -s)"     # start ssh-agent in bg
ssh-add ~/.ssh/id_rsa      # add your key to the ssh-agent

# Now copy the public key (~/.ssh/id_rsa.pub) exactly without adding newlines or whitespace to your 
clipboard!

# Again, on github.com, click you photo, goto settings / Add SSH-key, and add the key

ssh -T git@github.com # test your key
