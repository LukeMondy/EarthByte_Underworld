# EarthByte Additions - a toolbox for Underworld #

### What is this repository for? ###
The Underworld numerical modelling framework is a powerful tool for investigating geodynamic processes. On top of this, it is open source, which allows anyone to go in, view the code, modify it, and add to it. Furthermore, the Underworld team support this attitude by providing code that can generate 'toolboxes'. Toolboxes are essentially nice containers for additional code chunks to work within the Underworld framework.

The EarthByte_Additions toolbox is a collection of bits and pieces of code that have been developed within the EarthByte group at the University of Sydney School of Geosciences. Primarily it contains code that has been copied from the Underworld branch and modified in some way, though there are some chunks of original code. 

This work has been done by:

* Luke Mondy (PhD Candidate)

### How do I get set up? ###

* Install Underworld (probably from src)
* In the underworld src directory: 
```
#!bash

hg clone https://bitbucket.org/lmondy/earthbyte_additions
```

* Configure the toolbox
```
#!bash

UW_DIR=$PWD/../build ./configure.py <same config flags as you used for UW>
```

* Compile it
```
#!bash

UW_DIR=$PWD/../build ./scons.py
```

* In your Underworld input file, add: 

```
#!xml

    <list name="import" mergeType="merge">
        <param> earthbyte_additions </param>
    </list>
```


### Contribution guidelines ###

Feel free to submit bugs or new features.

### Who do I talk to? ###

* luke.s.mondy@gmail.com
* Underworld Developers
* You may also be interested in the Lithospheric Modelling Recipe