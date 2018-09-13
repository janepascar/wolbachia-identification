cd /N/dc2/scratch/chrichan/wolbachia/  

mkdir -p assembly_cov_med

cd assembly_cov_med

touch all_median_coverage.txt
rm all_median_coverage.txt
touch all_median_coverage.txt


for SAMPLE in Acromyrmex_echinatior Anopheles_gambiae Biorhiza_pallida_1 Biorhiza_pallida_2 Biorhiza_pallida_3 Callosobruchus_chinensis Ceratina_calcarata_1 Cynipini_1 Cynipini_2 Dactylopius_coccus Diabrotica_virgifera_virgifera_1 Diabrotica_virgifera_virgifera_2 Diabrotica_virgifera_virgifera_3 Diachasma_alloeum Diaphorina_citri_1 Diaphorina_citri_2 Diploeciton_nevermanni Diplolepis_spinosa_1 Diplolepis_spinosa_2 Drosophila_melanogaster_1 Drosophila_melanogaster_2 Drosophila_melanogaster_3 Drosophila_simulans_1 Drosophila_simulans_2 Drosophila_simulans_3 Drosophila_triauraria_1 Drosophila_yakuba_1 Drosophila_yakuba_2 Ecitophya_simulans Gerris_buenoi_1 Gerris_buenoi_2 Homalodisca_vitripennis_1 Homalodisca_vitripennis_2 Isocolus_centaureae_1 Isocolus_centaureae_2 Maconellicoccus_hirsutus Megacopta_cribraria Mycopsylla_fici_1 Mycopsylla_fici_2 Mycopsylla_proxima Pararge_aegeria Pediaspis_aceris_1 Pediaspis_aceris_2 Polygonia_c-album Pseudomyrmex_sp_PSW_54 Rhagoletis_zephyria Trichogramma_pretiosum Wuchereria_bancrofti Dirofilaria_immitis_1 Dirofilaria_immitis_2 Brugia_malayi Onchocerca_volvulus Onchocerca_ochengi Onchocerca_gutturosa Litomosoides_sigmodontis
do
    cat ../${SAMPLE}/round5/cov_data.txt | cut -f 2 | sort -n | awk -v sample=${SAMPLE} '{a[i++] = $1;} END {print sample "\t" a[int((i-1)/2)];}' >> all_median_coverage.txt
done

for SAMPLE in Bembidion_lapponicum Brugia_pahangi Delias_oraia Operophtera_brumata Rhagoletis_pomonella
do
    cat ../${SAMPLE}/round1/cov_data.txt | cut -f 2 | sort -n | awk -v sample=${SAMPLE} '{a[i++] = $1;} END {print sample "\t" a[int((i-1)/2)];}' >> all_median_coverage.txt
done

