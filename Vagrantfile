# -*- mode: ruby -*-
# vi: set ft=ruby :

# Vagrantfile API/syntax version. Don't touch unless you know what you're doing!
VAGRANTFILE_API_VERSION = "2"

Vagrant.configure(VAGRANTFILE_API_VERSION) do |config|
    # Set up the box
    config.vm.box = "ubuntu/xenial64"
    config.vm.provider "virtualbox" do |v|
      v.memory = 2048
      # v.cpus = 2
    end
    
    # Port forwarding
    # RStudio
    # Tunnel to localhost:8989 to avoid conflict with local RStudio server instance
    config.vm.network "forwarded_port", guest: 8787, host: 8888

    config.vm.hostname = "vagrant-dev-box"

    config.vm.provision :shell, path: "bootstrap.sh"

    # For R package development, uncomment these lines
    # replace /packageName with the R package name, 
    # then inside Rstudio create a new package in folder, 
    # move these around so that Rproj file is in /packageName
    #
    config.vm.synced_folder ".", "/vagrant", disabled: true
    config.vm.synced_folder ".", "/dremeR"


end
