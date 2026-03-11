process COLLECT_VERSIONS {

    label 'ultralow'

    input:
        path('versions_*')

    output:
        path("software_versions.yml"), emit: combined_versions

    script:
    """
    cat versions_* > all_versions.yml

    # Deduplicate: keep first occurrence of each tool entry
    python3 -c "
import sys

seen = set()
output = []
with open('all_versions.yml') as f:
    for line in f:
        stripped = line.strip()
        if stripped.endswith(':') and not stripped.startswith(' '):
            # Process header line - skip duplicates
            current_process = stripped
            if current_process not in seen:
                seen.add(current_process)
                output.append(line)
                is_dup = False
            else:
                is_dup = True
        elif stripped and not is_dup:
            output.append(line)

with open('software_versions.yml', 'w') as f:
    f.writelines(output)
" 2>/dev/null || cat all_versions.yml > software_versions.yml
    """
}
