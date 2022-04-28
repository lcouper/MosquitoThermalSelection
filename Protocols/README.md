# Introduction #

Using the Illumina DNA prep more efficiently (a 1/5 dilution of the protocol). 
Original protocol by Katie Solari. 
Slight formatting changes and explicit improvements for 96 well plate prep. See also Illumina protocol for reference. 

## Materials ##
 
- Genomic DNA (for yeast - I've been doing plate extractions using the Invitrogen Purelink 96 protocol)
- Illumina DNA prep
- Bead-linked transposase (BLT) - stored at 4°C 
- Tagmentation Buffer 1 (TB1) - stored at -20°C 
- Tagmentation Stop Buffer (TSB) - stored at 4°C 
- Tagmentation Wash Buffer (TWB) - stored at 4°C 
- Enhanced PCR Mix (EPM) - stored at -20°C 
- Sample Purification Beads (SPB) - stored at 4°C 
- Re-suspension Buffer (RSB) - stored at -20°C
- Illumina UDI primers (can also order identical primers from IDT)

** NOTE** make sure you keep track of which primers you use for each plate so that you can pool things appropriately for sequencing

## Procedure ##
DNA preparation stage 

1. Quantify DNA

2. Add water and gDNA to get at least 20ng of DNA in 6ul. More DNA is fine also, but all samples should be similar in concentration because this will impact the number of PCR cycles. [GRK note: I try to just keep samples within the ranges of the table below without trying to get them to be super close in concentration - I haven't verified how good this works yet]

Total DNA Input (ng) | 1-9 | 10-24 | 25-49 | 50-99 | 100-500 | 
--- | --- | --- | --- |--- |--- |--- |--- |--- |--- |--- |---
Number of PCR cycles | 12 | 8 | 6 | 5 | 5 | 

### Tagment Genomic DNA ###

3. Get Bead-linked Transposase (BLT) and Tagmentation Buffer 1 (TB1) to room temperature.

4. Vortex Bead-linked Transposases for 10 secondes to resuspend

5. Make Tagmentation Master mix:
216ul (2.25 x 96 samples) Bead-linked transposase
216ul (2.25 x 96 samples) Tagmentation Buffer 1
Vortex and Split into 12 PCR strip tubes (34ul per tube)

6. Add 4ul of Tagmentation Master Mix to each 6ul of gDNA and mix with pipette at least 10x.

7. Put in Thermocycler 
Lid 100°C
55°C for 15 minutes 
10°C hold

### Post Tagmentation Cleanup ###

8. Bring Tagmentation Stop Buffer and Tagmemtation Wash Buffer to room temperature

9. Add 20ul of Tagmentation Stop Buffer to each of 12 strip tubes

10. Add 2ul of Tagmentation Stop Buffer to each well of the sample plate, pipette up and down at least 10x to resuspend beads

11. Put in Thermocycler 
Lid 100°C
37°C for 15 minutes
10°C hold

12. Place sample plate on magnetic plate for 3 minutes and then remove and discard supernatent
NOTE: Easier to use multichannel here if you use the magnets that move beads to one side of the plate (rather than circular magnets)

13. Wash two times as follows:
a. Remove the sample plate form the magnet and add 20ul of Tagmentation Wash Buffer (TWB) to each well. Be sure to pipette up and down until beads are mixed.

NOTE: It can help to resuspend if you pipette down the wall that the beads are on a few times. Also, TWB is very soapy, so pipette slowly to avoid making bubbles. You can spin VERY briefly in centrifuge if needed to get rid of some bubbles. During this step, I also use 8 tips so it's easier to resuspend the beads (they are then all on the same wall).
b. Place on magnet for 3 minutes
c. Using multichannel and while the remove supernatent
NOTE: Don't worry too much about getting ALL the supernatent just yet...

14. Remove sample plate from the magnet, add 20ul of TWB to each well, resuspend beads

15. Seal and place back on magnet, keeping them in TWB until PCR amplification ready...

### Amplify Tagmented DNA ###
16. Thaw EPM and primer plate at room temperature. Invert to mix and then briefly centrifuge.

17. Make PCR Master Mix (pipette up and down to mix)
440ul (4.4ul x 100 samples) of EPM
440ul (4.4ul x 100 samples) of Nuclease Free water
 
18. Split PCR master mix into 8 strip tubes (109uL/ tube)

19. Go back to sample plate on the magnet,remove and discard supernatent while still on the magnet

20. Remove the sample plate from the magnet and add 8ul of PCR master mix, mix well with pipette and briefly centrifuge.

21. Vortex and briefly centrifuge primer plate(s).

22. Add 1ul adapter i7 and 1ul adapter i5 (or 2ul if using pre-mixed adapter plates from Illumina) to each tube. Vortex and briefly centrifuge.
NOTE: it's easier to remove the supernatent without disturbing the beads at this point if you use a P10 multichannel several times instead of P100. Here it's more important to get as much of the supernatent as you can.

23. Put in theromcycler

Lid 100°C
68°C for 3 minutes
98°C for 3 minutes
X cycles of: (see table at top of protocol for number of cycles depending on amount of gDNA you start with) 98°C for 45 seconds,
62°C for 30 seconds, 68°C for 2 min
68°C for 1 minute 10°C hold

### Clean Up Libraries ###
24. Bring Sample Purification Beads and Resuspension Buffer to room temperature. Vortex beads thoroughly to resuspend

25. Prepare fresh 80% ethanol

26. Remove sample plate from theromcycler. Centrifuge down and put on magnet for 5minutes.

27. Take 9 uL of supernatent and place in new plate

28. Take 4uL from each well and pool by column into strip tubes (8 samples pooled to make a 32uL pool for each column - you should now have a total of 12 strip tubes corresponding to each column of the plate).

29. Seal the plate with the remaining 5ul of each sample and save to return to later if necessary. Store at -20°C.

30. Add 28 uL of nuclease free water and 32 uL of Sample Purification Beads (SPB) to each 32uL pool and mix.

NOTE: This is a good stopping point. This is a DNA library before size selection!
NOTE: THE RATIO OF BEADS:(LIBRARY+WATER) is very important, as this sets the size of the fragments
 that are selected. If you end up pooling less than 8 samples together, be sure that you change the amount of water you add appropriately.

31. Sit at room temperature for 5 minutes and then place on the magnet for 5 minutes

32. Vortex Sample Purification Beads again and place 9uL of SPB into a new set of empty strip tubes

33. Add 75 uL of supernatent to new strip tube with 9uL SPB and mix thoroughly with pipette

34. Sit at room temperature for 5 minutes and then place on the magnet for 5 minutes

35. Remove and discard superntatent.

NOTE: Since we are doing some ethanol washes next, you don't need to worry too much about gettting every drop of the supernatent out.

36. Wash two times as follows:
a. Add 100ul 80% ethanol to each tube without mixing 
b. Let sit for 30 seconds
c. Without disturbing the beads, remove and discard supernatent

37. Remove residual ethanol with a P10 after second wash (you really want to try to get all you can here). 

38. Leave tubes open on magnet to air dry for exactly 5 minutes

39. Remove tubes from the magnet and add 22ul of Resuspension Buffer (RSB) to each tube. Resuspend with pipette.

NOTE: Leaving tubes air-drying for too long can cause the beads to dry out and crack which is bad

40. Sit at room temperature for 2 minutes then place on the magnet for 2 minutes 

41. While avoiding disturbing the beads, transfer 20ul of supernatent and put in new labeled tube - this is your end set of pooled libraries! Yay! Great job! (You can now quantify for further pooling, other downstream QC before submitting for sequencing).


