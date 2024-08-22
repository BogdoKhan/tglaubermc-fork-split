Open ROOT, on prompt write:

`root [0] .x 1_rootlogon.C`

After the compilation of shared libraries, the 

`root [1]`

should appear. You are ready to work; to proceed try following commands:

`root [1] runAndSaveNtuple(1000, "Xe", "W", 23.58)`

shall create .root file with the content of Ntuple, which you can check with TBrowser,
or:

```
root [1] TGlauberMC gmc ("Xe","W", 23.58)
root [2] gmc.Run(1000)
root [3] TNtuple* ntu = gmc.GetNtuple()
root [4] ntu->Draw("Npart")
```

will create this output

```
root [1] TGlauberMC gmc ("Xe","W", 23.58)
Setting up nucleus Xe
Setting up nucleus W
(TGlauberMC &) Name: Glauber_Xe_W Title: Glauber Xe+W Version
root [2] gmc.Run(1000)
Event # 900 x-sect = 5.6218 +- 0.187289 b
Done!
root [3] TNtuple* ntu = gmc.GetNtuple()
(TNtuple *) 0x7fffbedd18d0
root [4] ntu->Draw("Npart")
Info in <TCanvas::MakeDefCanvas>:  created default TCanvas with name c1
```
to try it with your own hands and draw specific ntuple in canvas.
