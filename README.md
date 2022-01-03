# t-hurdles
Code for "Do t-stat hurdles need to be raised?"

- Previous versions of code were in a zip file on my website
- New in this version:
  - Removed Fortran MEX stuff, that wasn't needed
  - Simplified / cleaned a lot of code
  - Uses the Chen-Zimmermann (CFR Forthcoming) open source dataset (instead of the data from the RAPS paper)
  - Allows for truncated student's t partial correlations to capture the fatter tails in correlations in the CFR paper
    - This richness doesn't change anything, but it's nice to check and pretty to look at
- Still to do: 
  - Finalize the robustness exhibit update