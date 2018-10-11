# FiniteSizeEffect
The codes here are for the paper "Finite size effects for spiking neural networks with spatially dependent coupling", which is posted on arXiv:
https://arxiv.org/pdf/1805.03601.pdf

The code now is ready for run on the newest version of Julia, v1.0.1. You can download at the following website:
https://julialang.org/downloads/

Source code is in the directory "src/". You can check the ".toml" files for the dependency versions. Using the module "FiniteSizeEffect", you can call all the convenient library in the old versions of Julia, and you won't need to worry about future update of Julia. 

To test the code, run "julia runcode.jl". Notice that for file "runcode.jl", you need to fill in your address of the file using push, otherwise julia will only load the default directory defined in "LOAD_PATH".
