def validateSampleSheetRows(rows, refCondition, design = 'condition') {
    if (!rows) {
        throw new IllegalArgumentException('The sample sheet contains no data rows.')
    }

    if (!(design in ['condition', 'donor_condition'])) {
        throw new IllegalArgumentException("Unsupported DESeq2 design '${design}'. Use 'condition' or 'donor_condition'.")
    }

    def requiredColumns = ['sample_id', 'fastq_1', 'fastq_2', 'condition']
    if (design == 'donor_condition') requiredColumns << 'donor'
    def missingColumns = requiredColumns.findAll { column -> !(column in rows[0].keySet()) }
    if (missingColumns) {
        throw new IllegalArgumentException("Sample sheet is missing required columns: ${missingColumns.join(', ')}")
    }

    def sampleIds = rows.collect { row -> row.sample_id?.trim() }
    if (rows.any { row -> row.sample_id != row.sample_id?.trim() }) {
        throw new IllegalArgumentException('Sample IDs must not have leading or trailing whitespace.')
    }
    if (sampleIds.any { sampleId -> !sampleId }) {
        throw new IllegalArgumentException('Sample IDs must not be empty.')
    }
    if (sampleIds.any { sampleId -> !(sampleId ==~ /[A-Za-z0-9][A-Za-z0-9_.-]*/) }) {
        throw new IllegalArgumentException('Sample IDs may contain only letters, numbers, dots, underscores, and hyphens, and must start with a letter or number.')
    }
    if (sampleIds.toSet().size() != sampleIds.size()) {
        throw new IllegalArgumentException('Sample IDs must be unique.')
    }

    def conditions = rows.collect { row -> row.condition?.trim() }
    if (rows.any { row -> row.condition != row.condition?.trim() }) {
        throw new IllegalArgumentException('Condition labels must not have leading or trailing whitespace.')
    }
    if (conditions.any { condition -> !condition }) {
        throw new IllegalArgumentException('Condition labels must not be empty.')
    }
    def conditionLevels = conditions.toSet()
    if (conditionLevels.size() != 2) {
        throw new IllegalArgumentException("Exactly two condition labels are required for the single DESeq2 contrast; found ${conditionLevels.size()}.")
    }
    if (!refCondition || !(refCondition in conditionLevels)) {
        throw new IllegalArgumentException("Reference condition '${refCondition}' is not present in the sample sheet: ${conditionLevels.join(', ')}.")
    }

    def samplesPerCondition = conditions.countBy { condition -> condition }
    def modelParameterCount = conditionLevels.size()
    if (design == 'donor_condition') {
        def donorCount = rows.collect { row -> row.donor?.trim() }.toSet().size()
        modelParameterCount = donorCount + 1
    }
    def residualDegreesOfFreedom = rows.size() - modelParameterCount
    if (residualDegreesOfFreedom < 1) {
        throw new IllegalArgumentException("The DESeq2 design '${design}' has no residual degrees of freedom (${rows.size()} samples, ${modelParameterCount} model parameters). Add biological replicates before analysis.")
    }
    def underReplicatedConditions = samplesPerCondition.findAll { _condition, count -> count < 2 }.keySet()
    if (underReplicatedConditions) {
        System.err.println("WARNING: Fewer than two samples are available for condition(s) ${underReplicatedConditions.join(', ')}. Replication is recommended for reliable dispersion estimation.")
    }

    if (design == 'donor_condition') {
        def donors = rows.collect { row -> row.donor?.trim() }
        if (rows.any { row -> row.donor != row.donor?.trim() }) {
            throw new IllegalArgumentException('Donor labels must not have leading or trailing whitespace.')
        }
        if (donors.any { donor -> !donor }) {
            throw new IllegalArgumentException('Donor labels must not be empty for the donor_condition design.')
        }
        if (donors.toSet().size() < 2) {
            throw new IllegalArgumentException('The donor_condition design requires at least two donors.')
        }
        def conditionsByDonor = [:].withDefault { [] as Set }
        rows.eachWithIndex { _row, index -> conditionsByDonor[donors[index]] << conditions[index] }
        def unpairedDonors = conditionsByDonor.findAll { _donor, labels -> labels != conditionLevels }.keySet()
        if (unpairedDonors) {
            throw new IllegalArgumentException("Every donor must have samples in both conditions for the donor_condition design; incomplete donors: ${unpairedDonors.join(', ')}.")
        }
    }

    rows.eachWithIndex { row, index ->
        ['fastq_1', 'fastq_2'].each { column ->
            def path = row[column]?.trim()
            if (!path) {
                throw new IllegalArgumentException("Sample '${sampleIds[index]}' has an empty ${column} path.")
            }
            if (!file(path).exists()) {
                throw new IllegalArgumentException("Sample '${sampleIds[index]}' ${column} file does not exist: ${path}")
            }
        }
    }

    return rows
}

def validateHisat2IndexFiles(indexFiles, indexPrefix) {
    if (!indexFiles) {
        throw new IllegalArgumentException("No HISAT2 index files found for prefix '${indexPrefix}'. Expected eight .ht2 or .ht2l files.")
    }

    def prefixName = indexPrefix.tokenize('/').last()
    def indexExtension = ['ht2', 'ht2l'].find { extension ->
        (1..8).every { part ->
            indexFiles.any { indexFile ->
                indexFile.fileName.toString() == "${prefixName}.${part}.${extension}"
            }
        }
    }

    if (!indexExtension) {
        throw new IllegalArgumentException("HISAT2 index '${indexPrefix}' is incomplete; expected ${prefixName}.1 through ${prefixName}.8 with a matching .ht2 or .ht2l suffix.")
    }

    return indexFiles.findAll { indexFile ->
        indexFile.fileName.toString().endsWith(".${indexExtension}")
    }
}
