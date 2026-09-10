#!/bin/bash
set -e

NGEO=3000
CDIR=${PWD}

for f in INCAR KPOINTS POTCAR; do
    if [[ ! -f "${CDIR}/${f}" ]]; then
        echo "ERROR: ${f} not found in ${CDIR}"
        exit 1
    fi
done

mkdir -p run

for i in $(seq 1 ${NGEO})
do
    ii=$(printf "%04d" ${i})

    if [[ ! -f "configs/POSCAR.${ii}" ]]; then
        echo "ERROR: configs/POSCAR.${ii} not found"
        exit 1
    fi

    mkdir -p "run/${ii}"

    cp "configs/POSCAR.${ii}" "run/${ii}/POSCAR"

    ln -sfn "${CDIR}/KPOINTS" "run/${ii}/KPOINTS"
    ln -sfn "${CDIR}/POTCAR"  "run/${ii}/POTCAR"

    if [[ ${i} -eq 1 ]]; then
        # First frame starts without a previous charge density
        cp "${CDIR}/INCAR" "run/${ii}/INCAR"

        if grep -Eq '^[[:space:]]*ICHARG[[:space:]]*=' "run/${ii}/INCAR"; then
            sed -i -E \
                's/^[[:space:]]*ICHARG[[:space:]]*=.*/ICHARG = 2/' \
                "run/${ii}/INCAR"
        else
            echo "ICHARG = 2" >> "run/${ii}/INCAR"
        fi
    else
        # Frames 0002-3000 use the normal INCAR (ICHARG = 1)
        ln -sfn "${CDIR}/INCAR" "run/${ii}/INCAR"
    fi

    if (( i % 100 == 0 || i == 1 || i == NGEO )); then
        echo "Initialized run/${ii}"
    fi
done

echo
echo "Initialization finished."
echo "run/0001/INCAR is a real file with ICHARG = 2."
echo "run/0002-run/3000 INCARs link to ${CDIR}/INCAR."