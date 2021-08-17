GenerateMomentMaps.py --> makes moments 0,1, & 2 as well as masked versions of these. 
CenterCubeOnCII.py --> reads in the combined cube and makes the central velocity 0 








# git pushes #
git add . 
git commit -m "date"
git push | git push -u origin master

# Accidently push a large file # 
git reset --soft HEAD~1 # to roll back a commit
add files to the .gitignore file 

# Adding the Oath Token. 
Git token from github under developer settings
git remote remove origin 
git remote add origin https://<token>@github.com/TrystanScottLambert/Hz7_ISM.git
