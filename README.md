### AMBER-PMEMD builder for Sunway TaihuLight

We discarded the original build system of Amber, instead, we migrated a part of building system of CESM for better testing support.

pmemd_min is the minimum building system for PMEMD, use ./create_case.py to create new test cases.

### This repo also contains optimimzed PMEMD code for Sunway TaihuLight's SW26010 processor.
SourceMods contains the newest version of our code, due to Amber is not an open source software, we only published our added files.

In addition, for our modified files, corresponding object files is also in SourceMods.

