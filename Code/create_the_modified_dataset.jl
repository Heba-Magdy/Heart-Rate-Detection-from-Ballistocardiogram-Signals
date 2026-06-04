function copy_modified_dataset()

    src_root = raw"C:\Users\Heba\Desktop\PhD_Courses\AS\project\dataset\dataset\data"
    dst_root = raw"C:\Users\Heba\Desktop\PhD_Courses\AS\project\dataset\dataset\modified_signals"

    mkpath(dst_root)

    patients = readdir(src_root)

    for patient in patients

        patient_path = joinpath(src_root, patient)

        bcg_path = joinpath(patient_path, "BCG")
        rr_path  = joinpath(patient_path, "Reference", "RR")

        if !isdir(bcg_path) || !isdir(rr_path)
            println("Skipping $patient (missing folders)")
            continue
        end

        bcg_files = readdir(bcg_path)
        rr_files  = readdir(rr_path)

        if isempty(bcg_files) || isempty(rr_files)
            println("Skipping $patient (empty folders)")
            continue
        end

        copied_any = false

        # ==========================================
        # LOOP OVER ALL BCG FILES
        # ==========================================
        for bcg_file in bcg_files

            bcg_base = splitext(bcg_file)[1]
            bcg_base = bcg_base[1:length(bcg_base)-3]
            rr_match = ""

            # find matching RR file
            for f in rr_files
                rr_base = splitext(f)[1]

                if occursin(bcg_base, rr_base) || occursin(rr_base, bcg_base)
                    rr_match = f
                    break
                end
            end

            # skip if no match
            if rr_match == ""
                println("No RR match for $patient → $bcg_file")
                continue
            end

            # --------------------------
            # COPY BOTH FILES
            # --------------------------
            src_bcg = joinpath(bcg_path, bcg_file)
            src_rr  = joinpath(rr_path, rr_match)

            dst_patient = joinpath(dst_root, patient)
            dst_bcg_dir = joinpath(dst_patient, "BCG")
            dst_rr_dir  = joinpath(dst_patient, "RR")

            mkpath(dst_bcg_dir)
            mkpath(dst_rr_dir)

            cp(src_bcg, joinpath(dst_bcg_dir, bcg_file); force=true)
            cp(src_rr,  joinpath(dst_rr_dir, rr_match); force=true)

            println("Copied $patient → $bcg_file ↔ $rr_match")

            copied_any = true
        end

        if !copied_any
            println("Skipping $patient (no valid pairs found)")
        end
    end
end

copy_modified_dataset()