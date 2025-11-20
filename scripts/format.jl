using JuliaFormatter

"""
    format_all_julia_files(directory::String)

Formate récursivement tous les fichiers Julia (.jl) dans le dossier spécifié.
"""
function format_all_julia_files(directory::String)
    # Parcourir récursivement le dossier
    for (root, dirs, files) in walkdir(directory)
        for file in files
            if endswith(file, ".jl")
                filepath = joinpath(root, file)
                println("Formatting $filepath")
                format_file(filepath; overwrite=true)
            end
        end
    end
end
