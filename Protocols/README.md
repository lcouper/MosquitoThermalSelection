# Introduction #

Using the Illumina DNA prep more efficiently (a 1/5 dilution of the protocol). 
Original protocol by Katie Solari. 
Slight formatting changes and explicit improvements for 96 well plate prep by GRK. See also Illumina protocol for reference. Materials
 
• Genomic DNA (for yeast - I've been doing plate extractions using the Invitrogen Purelink 96 protocol)
• Illumina DNA prep
oBead-linked transposase (BLT) - stored at 4°C oTagmentation Buffer 1 (TB1) - stored at -20°C oTagmentation Stop Buffer (TSB) - stored at 4°C oTagmentation Wash Buffer (TWB) - stored at 4°C oEnhanced PCR Mix (EPM) - stored at -20°C oSample Purification Beads (SPB) - stored at 4°C oRe-suspension Buffer (RSB) - stored at -20°C
• Illumina UDI primers (can also order identical primers from IDT)

**NOTE** make sure you keep track of which primers you use for each plate so that you can pool things appropriately for sequencing

Procedure
DNA preparation stage 
1. Quantify DNA
2. Add water and gDNA to get at least 20ng of DNA in 6ul. More DNA is fine also, but all samples should be similar in concentration because this will impact the number of PCR cycles. [GRK note: I try to just keep samples within the ranges of the table below without trying to get them to be super close in concentration - I haven't verified how good this works yet]
1 2 3 4 5 6
C
3.
A
B
    
Total DNA Input (ng)
Number of PCR cycles
    
1-9
12
     
 
10-24
8
      
25-49
6
     
 
50-99
5
     
 
100-500
5
 Number of PCR cycles based on Total DNA Input
 Tagment Genomic DNA
3. Get Bead-linked Transposase (BLT) and Tagmentation Buffer 1 (TB1) to room temperature.
4. Vortex Bead-linked Transposases for 10 secondes to resuspend
5. Make Tagmentation Master mix:
•
• • •
6. Add 4ul of Tagmentation Master Mix to each 6ul of gDNA and mix with pipette at least 10x.
7. Put in Thermocycler (saved under GRK/ILLUM_1)
 • • • •
Lid 100°C
55°C for 15 minutes 10°C hold
Post Tagmentation Cleanup
216ul (2.25 x 96 samples) Bead-linked transposase 216ul (2.25 x 96 samples) Tagmentation Buffer 1 Vortex and Split into 12 PCR strip tubes (34ul per tube)
 8. Bring Tagmentation Stop Buffer and Tagmemtation Wash Buffer to room temperature
9. Add 20ul of Tagmentation Stop Buffer to each of 12 strip tubes
10. Add2ulofTagmentationStopBuffertoeachwellofthesampleplate,pipetteupanddownatleast10xto
 resuspend beads
11. Put in Thermocycler (saved under GRK/ILLUM_2)
• •
Lid 100°C
37°C for 15 minutes
• •
13. Wash two times as follows:
•
 10°C hold
12. Place sample plate on magnetic plate for 3 minutes and then remove and discard supernatent
NOTE: Easier to use multichannel here if you use the magnets that move beads to one side of the plate (rather than circular magnets)
a. Remove the sample plate form the magnet and add 20ul of Tagmentation Wash Buffer (TWB) to each well. Be sure to pipette up and down until beads are mixed.

•
• • •
NOTE: It can help to resuspend if you pipette down the wall that the beads are on a few times. Also, TWB is very soapy, so pipette slowly to avoid making bubbles. You can spin VERY briefly in centrifuge if needed to get rid of some bubbles. During this step, I also use 8 tips so it's easier to resuspend the beads (they are then all on the same wall).
 b. Place on magnet for 3 minutes
c. Using multichannel and while the remove supernatent
NOTE: Don't worry too much about getting ALL the supernatent just yet...
 14. Remove sample plate from the magnet, add 20ul of TWB to each well, resuspend beads
15. Seal and place back on magnet, keeping them in TWB until PCR amplification ready...
Amplify Tagmented DNA
16. ThawEPMandprimerplateonice(Ifindthatitneverthawsonice,soIjustthawitatroomtemperature). Invert to mix and then briefly centrifuge.
17. MakePCRmastermix(pipetteupanddowntomix):
•
•
•
18. SplitPCRmastermixinto8striptubes(109ul/tube)
440ul (4.4ul x 100 samples) of EPM
440ul (4.4ul x 100 samples) of Nuclease Free water
 19. Gobacktosampleplateonthemagnet,removeanddiscardsupernatentwhilestillonthemagnet.
•
20. Remove the sample plate from the magnet and add 8ul of PCR master mix, mix well with pipette and briefly centrifuge.
21. Vortex and briefly centrifuge primer plate(s).
22. Add 1ul adapter i7 and 1ul adapter i5 (or 2ul if using pre-mixed adapter plates from Illumina) to each tube. Vortex and briefly centrifuge.
NOTE: it's easier to remove the supernatent without disturbing the beads at this point if you use a P10 multichannel several times instead of P100. Here it's more important to get as much of the supernatent as you can.
 23. Putintheromcycler(savedunderGRK/ILLUM_3)
• •
Lid 100°C
68°C for 3 minutes
• • •
98°C for 3 minutes
X cycles of: (see table at top of protocol for number of cycles depending on amount of gDNA you start with) 98°C for 45 seconds,
 
• • • • •
62°C for 30 seconds, 68°C for 2 min
68°C for 1 minute 10°C hold
Clean Up Libraries
  24. BringSamplePurificationBeadsandResuspensionBufffertoroomtemperature.Vortexbeadsthoroughlyto resuspend.
25. Preparefresh80%ethanol
26. Removesampleplatefromtheromcycler.Centrifugedownandputonmagnetfor5minutes.
27. Take9ulofsupernatentandplaceinnewplate.
•
28. Take4ulfromeachwellandpoolbycolumnintostriptubes(8samplespooledtomakea32ulpoolforeach column - you should now have a total of 12 strip tubes corresponding to each column of the plate).
29. Seal the plate with the remaining 5ul of each sample and save to return to later if necessary. Store at -20°C.
30. Add28ulofnucleasefreewaterand32ulofSamplePurificationBeads(SPB)toeach32ulpoolandmix.
•
NOTE: This is a good stopping point. This is a DNA library before size selection!
NOTE: THE RATIO OF BEADS:(LIBRARY+WATER) is very important, as this sets the size of the fragments
 that are selected. If you end up pooling less than 8 samples together, be sure that you change the amount of water you add appropriately.
31. Sitatroomtemperaturefor5minutesandthenplaceonthemagnetfor5minutes.
32. VortexSamplePurificaitonBeadsagainandplace9ulofSPBintoanewsetofemptystriptubes.
33. Add75ulofsupernatenttonewstriptubewith9ulSPBandmixthoroughlywithpipette.
•
NOTE: Katie's original protocol calls for 25ul supernatent and 3ul SPB here. This results in slightly lower yield, so I changed it to keep more of the DNA. Again, ratio is super important here.
 34. Sitatroomtemperaturefor5minutesandthenplaceonthemagnetfor5minutes. 35. Remove and discard superntatent.
•
NOTE: Since we are doing some ethanol washes next, you don't need to worry too much about gettting every drop of the supernatent out.
 36. Wash two times as follows:
• •
a. Add 100ul 80% ethanol to each tube without mixing b. Let sit for 30 seconds

•
•
39. Remove tubes from the magnet and add 22ul of Resuspension Buffer (RSB) to each tube. Resuspend with pipette.
c. Without disturbing the beads, remove and discard supernatent
37. Remove residual ethanol with a P10 after second wash (you really want to try to get all you can here). 38. Leave tubes open on magnet to air dry for exactly 5 minutes
NOTE: Leaving tubes air-drying for too long can cause the beads to dry out and crack which is bad
40. Sitatroomtemperaturefor2minutesandthenplaceonthemagnetfor2minutes.
41. While avoiding disturbing the beads, transfer 20ul of supernatent and put in new labeled tube - this is your end set of pooled libraries! Yay! Great job! (You can now quantify for further pooling, other downstream QC before submitting for sequencing).
•
NOTE: Katie's original protocol originally calls for resuspending in 12ul RSB and getting 12ul of supernatent. I found it was quite hard to avoid getting any beads in the pipette, so I just added more RSB and abandon a couple of microliters to avoid the headache of struggling with getting beads. If you do end up pipetting some beads, pipette it back into the tube, wait a few minutes, and try again.
