          seed =  1

       seqfile = /global/scratch2/rohitkolora/Rockfish/Genomes/trees/nucleartree/BPP/ALL/rockfish.phy.txt.mod
      Imapfile = /global/scratch2/rohitkolora/Rockfish/Genomes/trees/nucleartree/BPP/ALL/rockfish.lmod.txt
       outfile = /global/scratch2/rohitkolora/Rockfish/Genomes/trees/nucleartree/BPP/ALL/A01/r3/out3.txt
      mcmcfile = /global/scratch2/rohitkolora/Rockfish/Genomes/trees/nucleartree/BPP/ALL/A01/r3/mcmc3.txt

  speciesdelimitation = 0 * fixed species tree
* speciesdelimitation = 1 0 2    * species delimitation rjMCMC algorithm0 and finetune(e)
* speciesdelimitation = 1 1 2 1 * species delimitation rjMCMC algorithm1 finetune (a m)
         speciestree = 1 * 0.4 0.2 0.1   * speciestree pSlider ExpandRatio ShrinkRatio

   speciesmodelprior = 1  * 0: uniform LH; 1:uniform rooted trees; 2: uniformSLH; 3: uniformSRooted

*  species&tree 
  species&tree = 83  Sebastolobus_alascanus Sebastolobus_altivelis Sebastiscus_albofasciatus Sebastiscus_tertius Adelosebastes_latens Hozukius_emblemarius Hozukius_guyotensis Sebastes_aleutianus Sebastes_baramenuke Sebastes_alutus Sebastes_ciliatus Sebastes_variabilis Sebastes_polyspinis Sebastes_crameri Sebastes_reedi Sebastes_auriculatus Sebastes_dalli Sebastes_rastrelliger Sebastes_aurora Sebastes_caurinus Sebastes_carnatus Sebastes_maliger Sebastes_atrovirens Sebastes_nebulosus Sebastes_saxicola Sebastes_semicinctus Sebastes_elongatus Sebastes_chlorostictus Sebastes_rosenblatti Sebastes_constellatus Sebastes_oculatus Sebastes_exsul Sebastes_helvomaculatus Sebastes_umbrosus Sebastes_ensifer Sebastes_rosaceus Sebastes_melanostomus Sebastes_miniatus Sebastes_pinniger Sebastes_babcocki Sebastes_diploproa Sebastes_rubrivinctus Sebastes_serriceps Sebastes_nigrocinctus Sebastes_levis Sebastes_goodei Sebastes_jordani Sebastes_paucispinis Sebastes_ruberrimus Sebastes_hopkinsi Sebastes_moseri Sebastes_proriger Sebastes_variegatus Sebastes_wilsoni Sebastes_zacentrus Sebastes_diaconus Sebastes_mystinus Sebastes_entomelas Sebastes_flavidus Sebastes_melanops Sebastes_fasciatus Sebastes_mentella Sebastes_glaucus Sebastes_itinus Sebastes_minor Sebastes_steindachneri Sebastes_hubbsi Sebastes_oblongus Sebastes_koreanus Sebastes_nudus Sebastes_pachycephalus Sebastes_nivosus Sebastes_schlegelii Sebastes_taczanowskii Sebastes_trivittatus Sebastes_zonatus Sebastes_inermis Sebastes_joyneri Sebastes_thompsoni Sebastes_scythropus Sebastes_kiyomatsui Sebastes_iracundus Sebastes_matsubarae
                   1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
                 ((Sebastolobus_alascanus, Sebastolobus_altivelis), ((Sebastiscus_albofasciatus, Sebastiscus_tertius), (Adelosebastes_latens, ((Hozukius_emblemarius, Hozukius_guyotensis), (((((((Sebastes_aleutianus, Sebastes_baramenuke), (Sebastes_alutus, ((Sebastes_ciliatus, Sebastes_variabilis), Sebastes_polyspinis))), (Sebastes_crameri, Sebastes_reedi)), (((((((((((Sebastes_auriculatus, (Sebastes_dalli, Sebastes_rastrelliger)), ((((Sebastes_aurora, Sebastes_caurinus), (Sebastes_carnatus, Sebastes_maliger)), Sebastes_atrovirens), Sebastes_nebulosus)), (Sebastes_saxicola, Sebastes_semicinctus)), Sebastes_elongatus), ((((Sebastes_chlorostictus, Sebastes_rosenblatti), (((Sebastes_constellatus, Sebastes_oculatus), Sebastes_exsul), (Sebastes_helvomaculatus, Sebastes_umbrosus))), Sebastes_ensifer), Sebastes_rosaceus)), Sebastes_melanostomus), (Sebastes_miniatus, Sebastes_pinniger)), ((((Sebastes_babcocki, (Sebastes_diploproa, ((Sebastes_rubrivinctus, Sebastes_serriceps), Sebastes_nigrocinctus))), Sebastes_levis), (Sebastes_goodei, (Sebastes_jordani, Sebastes_paucispinis))), Sebastes_ruberrimus)), (Sebastes_hopkinsi, Sebastes_moseri)), (Sebastes_proriger, (Sebastes_variegatus, (Sebastes_wilsoni, Sebastes_zacentrus)))), (((Sebastes_diaconus, Sebastes_mystinus), Sebastes_entomelas), (Sebastes_flavidus, Sebastes_melanops)))), ((((Sebastes_fasciatus, Sebastes_mentella), Sebastes_glaucus), (Sebastes_itinus, (Sebastes_minor, Sebastes_steindachneri))), ((((((Sebastes_hubbsi, Sebastes_oblongus), (Sebastes_koreanus, (Sebastes_nudus, Sebastes_pachycephalus))), Sebastes_nivosus), ((Sebastes_schlegelii, Sebastes_taczanowskii), (Sebastes_trivittatus, Sebastes_zonatus))), ((Sebastes_inermis, Sebastes_joyneri), Sebastes_thompsoni)), (Sebastes_scythropus, Sebastes_kiyomatsui)))), Sebastes_iracundus), Sebastes_matsubarae)))));
         diploid =   1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
   
       usedata = 1  * 0: no data (prior); 1:seq like
         nloci = 624 * number of data sets in seqfile

     cleandata = 1    * remove sites with ambiguity data (1:yes, 0:no)?

    thetaprior = 22 0.05 E  # invgamma(a, b) for theta
      tauprior = 3 0.14   # invgamma(a, b) for root tau & Dirichlet(a) for other tau's

*     heredity = 1 4 4
*    locusrate = 1 5

*      finetune =  1: 5 0.001 0.001  0.001 0.3 0.33 1.0  # finetune for GBtj, GBspr, theta, tau, mix, locusrate, seqerr

         print = 1 0 0 0   * MCMC samples, locusrate, heredityscalars, Genetrees
        burnin = 1000
      sampfreq = 2
       nsample = 150000
       threads = 32

