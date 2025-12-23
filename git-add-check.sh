#!/bin/bash
#
# Git add wrapper that checks for files >100MB and excludes them
#
# INSTALLATION:
# Add the following line to your ~/.bashrc or ~/.zshrc:
#
#   source /home/afreiburger/repos/codiffusion_bioreactor/git-add-check.sh
#
# This will override 'git add' to automatically check file sizes.
#

git() {
    if [ "$1" = "add" ]; then
        shift  # Remove 'add' from arguments
        __git_add_with_size_check "$@"
    else
        command git "$@"
    fi
}

__git_add_with_size_check() {
    # Maximum file size in bytes (100MB)
    local MAX_SIZE=$((100 * 1024 * 1024))
    local MAX_SIZE_HUMAN="100MB"

    # Arrays to track files
    local FILES_TO_ADD=()
    local LARGE_FILES=()
    local LARGE_FILE_DETAILS=()
    local GIT_ARGS=()
    local FILE_ARGS=()

    # Separate git options from file arguments
    for arg in "$@"; do
        if [[ "$arg" == -* ]]; then
            GIT_ARGS+=("$arg")
        else
            FILE_ARGS+=("$arg")
        fi
    done

    # If no file args, use current directory
    if [ ${#FILE_ARGS[@]} -eq 0 ]; then
        FILE_ARGS=(".")
    fi

    # Function to check a single file
    __check_file() {
        local file="$1"
        local FILE_SIZE

        if [ -f "$file" ]; then
            FILE_SIZE=$(stat -c%s "$file" 2>/dev/null || stat -f%z "$file" 2>/dev/null)
            if [ "$FILE_SIZE" -gt "$MAX_SIZE" ]; then
                local FILE_SIZE_MB=$(echo "scale=2; $FILE_SIZE / 1024 / 1024" | bc)
                LARGE_FILES+=("$file")
                LARGE_FILE_DETAILS+=("$file (${FILE_SIZE_MB}MB)")
                return 1
            fi
        fi
        return 0
    }

    # Process file arguments
    for pattern in "${FILE_ARGS[@]}"; do
        if [ -d "$pattern" ]; then
            # It's a directory - find all files
            while IFS= read -r -d '' file; do
                if __check_file "$file"; then
                    FILES_TO_ADD+=("$file")
                fi
            done < <(find "$pattern" -type f -print0 2>/dev/null)
        elif [ -f "$pattern" ]; then
            # It's a single file
            if __check_file "$pattern"; then
                FILES_TO_ADD+=("$pattern")
            fi
        else
            # Could be a glob pattern or doesn't exist - let git handle it
            # But first try to expand it
            local expanded=($pattern)
            for item in "${expanded[@]}"; do
                if [ -f "$item" ]; then
                    if __check_file "$item"; then
                        FILES_TO_ADD+=("$item")
                    fi
                elif [ -e "$item" ]; then
                    FILES_TO_ADD+=("$item")
                fi
            done
            # If nothing expanded, pass it to git as-is
            if [ ${#expanded[@]} -eq 0 ] || [ ! -e "${expanded[0]}" ]; then
                FILES_TO_ADD+=("$pattern")
            fi
        fi
    done

    # Report large files if found
    if [ ${#LARGE_FILES[@]} -gt 0 ]; then
        echo ""
        echo "=========================================="
        echo "  WARNING: Large files EXCLUDED (>$MAX_SIZE_HUMAN)"
        echo "=========================================="
        echo ""
        echo "The following files were NOT added because they exceed $MAX_SIZE_HUMAN:"
        echo ""
        for detail in "${LARGE_FILE_DETAILS[@]}"; do
            echo "  - $detail"
        done
        echo ""
        echo "File details:"
        for file in "${LARGE_FILES[@]}"; do
            ls -lh "$file" 2>/dev/null
        done
        echo ""
        echo "Consider:"
        echo "  - Adding to .gitignore: echo '${LARGE_FILES[0]}' >> .gitignore"
        echo "  - Using Git LFS: git lfs track '*.${LARGE_FILES[0]##*.}'"
        echo "=========================================="
        echo ""
    fi

    # Add the remaining files
    if [ ${#FILES_TO_ADD[@]} -gt 0 ]; then
        command git add "${GIT_ARGS[@]}" "${FILES_TO_ADD[@]}"
        local result=$?
        if [ $result -eq 0 ]; then
            echo "Successfully staged ${#FILES_TO_ADD[@]} file(s)."
            if [ ${#LARGE_FILES[@]} -gt 0 ]; then
                echo "(${#LARGE_FILES[@]} large file(s) excluded)"
            fi
        fi
        return $result
    else
        if [ ${#LARGE_FILES[@]} -gt 0 ]; then
            echo "No files staged - all matched files exceed the size limit."
        else
            echo "No files to add."
        fi
        return 0
    fi
}

echo "Git add wrapper loaded: files >100MB will be automatically excluded."
