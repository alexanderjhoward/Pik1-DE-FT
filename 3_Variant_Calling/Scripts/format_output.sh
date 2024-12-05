#!/bin/bash
# Collect formatted sequences
cat Variants/N22/*_E1.pep | awk '/^>/ { print (NR==1 ? "" : RS) $0; next } { printf "%s", $0 } END { printf RS }' > Output/N22/N22_E1.pep
cat Variants/N22/*_E2.pep | awk '/^>/ { print (NR==1 ? "" : RS) $0; next } { printf "%s", $0 } END { printf RS }' > Output/N22/N22_E2.pep
cat Variants/N22/*_E1.fa | awk '/^>/ { print (NR==1 ? "" : RS) $0; next } { printf "%s", $0 } END { printf RS }' > Output/N22/N22_E1.fa
cat Variants/N22/*_E2.fa | awk '/^>/ { print (NR==1 ? "" : RS) $0; next } { printf "%s", $0 } END { printf RS }' > Output/N22/N22_E2.fa
cat Variants/Vars/*_P1E1.pep | awk '/^>/ { print (NR==1 ? "" : RS) $0; next } { printf "%s", $0 } END { printf RS }' > Output/Vars/Vars_P1E1.pep
cat Variants/Vars/*_P1E2.pep | awk '/^>/ { print (NR==1 ? "" : RS) $0; next } { printf "%s", $0 } END { printf RS }' > Output/Vars/Vars_P1E2.pep
cat Variants/Vars/*_P1E3.pep | awk '/^>/ { print (NR==1 ? "" : RS) $0; next } { printf "%s", $0 } END { printf RS }' > Output/Vars/Vars_P1E3.pep
cat Variants/Vars/*_P2E1.pep | awk '/^>/ { print (NR==1 ? "" : RS) $0; next } { printf "%s", $0 } END { printf RS }' > Output/Vars/Vars_P2E1.pep
cat Variants/Vars/*_P2E2.pep | awk '/^>/ { print (NR==1 ? "" : RS) $0; next } { printf "%s", $0 } END { printf RS }' > Output/Vars/Vars_P2E2.pep
cat Variants/Vars/*_P1E1.fa | awk '/^>/ { print (NR==1 ? "" : RS) $0; next } { printf "%s", $0 } END { printf RS }' > Output/Vars/Vars_P1E1.fa
cat Variants/Vars/*_P1E2.fa | awk '/^>/ { print (NR==1 ? "" : RS) $0; next } { printf "%s", $0 } END { printf RS }' > Output/Vars/Vars_P1E2.fa
cat Variants/Vars/*_P1E3.fa | awk '/^>/ { print (NR==1 ? "" : RS) $0; next } { printf "%s", $0 } END { printf RS }' > Output/Vars/Vars_P1E3.fa
cat Variants/Vars/*_P2E1.fa | awk '/^>/ { print (NR==1 ? "" : RS) $0; next } { printf "%s", $0 } END { printf RS }' > Output/Vars/Vars_P2E1.fa
cat Variants/Vars/*_P2E2.fa | awk '/^>/ { print (NR==1 ? "" : RS) $0; next } { printf "%s", $0 } END { printf RS }' > Output/Vars/Vars_P2E2.fa

# Format to CSV
cat Output/N22_E1.pep | sed 's/>//' | paste -d ","  - - | sed '1s/^/id,seq\n/' | sed 's/_1//g' | awk '{ sub("\r$", ""); print }' > Output/N22/N22_E1_pep.csv
cat Output/N22_E2.pep | sed 's/>//' | paste -d ","  - - | sed '1s/^/id,seq\n/' | sed 's/_1//g' | awk '{ sub("\r$", ""); print }' > Output/N22/N22_E2_pep.csv
cat Output/N22_E1.fa | sed 's/>//' | paste -d ","  - - | sed '1s/^/id,seq\n/' | sed 's/_1//g' | awk '{ sub("\r$", ""); print }' > Output/N22/N22_E1_dna.csv
cat Output/N22_E2.fa | sed 's/>//' | paste -d ","  - - | sed '1s/^/id,seq\n/' | sed 's/_1//g' | awk '{ sub("\r$", ""); print }' > Output/N22/N22_E2_dna.csv
cat Output/Vars_P1E1.pep | sed 's/>//' | paste -d ","  - - | sed '1s/^/id,seq\n/' | sed 's/_1//g' | awk '{ sub("\r$", ""); print }' > Output/Vars/Vars_P1E1_pep.csv
cat Output/Vars_P1E2.pep | sed 's/>//' | paste -d ","  - - | sed '1s/^/id,seq\n/' | sed 's/_1//g' | awk '{ sub("\r$", ""); print }' > Output/Vars/Vars_P1E2_pep.csv
cat Output/Vars_P1E3.pep | sed 's/>//' | paste -d ","  - - | sed '1s/^/id,seq\n/' | sed 's/_1//g' | awk '{ sub("\r$", ""); print }' > Output/Vars/Vars_P1E3_pep.csv
cat Output/Vars_P2E1.pep | sed 's/>//' | paste -d ","  - - | sed '1s/^/id,seq\n/' | sed 's/_1//g' | awk '{ sub("\r$", ""); print }' > Output/Vars/Vars_P2E1_pep.csv
cat Output/Vars_P2E2.pep | sed 's/>//' | paste -d ","  - - | sed '1s/^/id,seq\n/' | sed 's/_1//g' | awk '{ sub("\r$", ""); print }' > Output/Vars/Vars_P2E2_pep.csv
cat Output/Vars_P1E1.fa | sed 's/>//' | paste -d ","  - - | sed '1s/^/id,seq\n/' | sed 's/_1//g' | awk '{ sub("\r$", ""); print }' > Output/Vars/Vars_P1E1_dna.csv
cat Output/Vars_P1E2.fa | sed 's/>//' | paste -d ","  - - | sed '1s/^/id,seq\n/' | sed 's/_1//g' | awk '{ sub("\r$", ""); print }' > Output/Vars/Vars_P1E2_dna.csv
cat Output/Vars_P1E3.fa | sed 's/>//' | paste -d ","  - - | sed '1s/^/id,seq\n/' | sed 's/_1//g' | awk '{ sub("\r$", ""); print }' > Output/Vars/Vars_P1E3_dna.csv
cat Output/Vars_P2E1.fa | sed 's/>//' | paste -d ","  - - | sed '1s/^/id,seq\n/' | sed 's/_1//g' | awk '{ sub("\r$", ""); print }' > Output/Vars/Vars_P2E1_dna.csv
cat Output/Vars_P2E2.fa | sed 's/>//' | paste -d ","  - - | sed '1s/^/id,seq\n/' | sed 's/_1//g' | awk '{ sub("\r$", ""); print }' > Output/Vars/Vars_P2E2_dna.csv
