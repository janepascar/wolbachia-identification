cd /N/dc2/scratch/chrichan/wolbachia

module load r

for SPECIES in Acromyrmex_echinatior Anopheles_gambiae Bembidion_lapponicum Biorhiza_pallida_1 Biorhiza_pallida_2 Biorhiza_pallida_3 Callosobruchus_chinensis Ceratina_calcarata_1 Cynipini_1 Cynipini_2 Dactylopius_coccus Delias_oraia Diabrotica_virgifera_virgifera_1 Diabrotica_virgifera_virgifera_2 Diabrotica_virgifera_virgifera_3 Diachasma_alloeum Diaphorina_citri_1 Diaphorina_citri_2 Diploeciton_nevermanni Diplolepis_spinosa_1 Diplolepis_spinosa_2 Drosophila_melanogaster_1 Drosophila_melanogaster_2 Drosophila_melanogaster_3 Drosophila_simulans_1 Drosophila_simulans_2 Drosophila_simulans_3 Drosophila_triauraria_1 Drosophila_yakuba_1 Drosophila_yakuba_2 Ecitophya_simulans Gerris_buenoi_1 Gerris_buenoi_2 Homalodisca_vitripennis_1 Homalodisca_vitripennis_2 Isocolus_centaureae_1 Isocolus_centaureae_2 Maconellicoccus_hirsutus Megacopta_cribraria Mycopsylla_fici_1 Mycopsylla_fici_2 Mycopsylla_proxima Pararge_aegeria Pediaspis_aceris_1 Pediaspis_aceris_2 Polygonia_c-album Pseudomyrmex_sp_PSW_54 Rhagoletis_pomonella Rhagoletis_zephyria Trichogramma_pretiosum
do
    cd ${SPECIES}
    Rscript --vanilla ../plot_histo.R
    cp coverage_histogram.pdf ../pdfs/${SPECIES}.pdf
    cd ..
done

#For species that have very high variance in Wolbachia coverage, plotting
# with the x-axis on a log scale produces a nicer histogram
for SPECIES in Operophtera_brumata
do
    cd ${SPECIES}
    Rscript --vanilla ../plot_histo_log.R
    cp coverage_histogram.pdf ../pdfs/${SPECIES}.pdf
    cd ..
done
