message("pdir is required to already be on github")

#  git clone https://github.com/Jadamso/GeoCleanR.git

#-------------------------------------------------------------------
##################
# Upload Package
################## 

#"git add . && git ci && git push"
repo <- git2r::repository(pdir)
git2r::add(repo, ".")
git2r::commit(repo, "new R commit")

#push(repo, credentials= cred_token(token = "GITHUB_PAT"))
system( paste0("cd ", pdir, "&& git push"))


#-------------------------------------------------------------------
##################
# Install Package Locally 
################## 

devtools::install(packg) ## Works

#-------------------------------------------------------------------
##################
# Install Package From Github
################## 
## Public Package

#install.packages("devtools")
devtools::install_github( paste0("Jadamso/",packg), subdir=packg)
citation(pack)
    

print("Done")


