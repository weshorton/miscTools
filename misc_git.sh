###
### Do a sparse checkout of a single sub-directory within a repo
###

# Initialize repo
git init <repo>; cd <repo>

# Add Origin
git remote add origin <url>

# Enable sparse checkout in config
git config core.sparsecheckout true

# Add desired directories to sparse-checkout file
echo "dirName/" >> .git/info/sparse-checkout

# Option 1. Pull from master
git pull origin master

# Option 2. Use shallow clone to cut off history (instead of pulling entire commit history)
git pull --depth=1 origin master
