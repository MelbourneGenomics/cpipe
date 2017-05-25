from email.utils import parseaddr
from pandas_schema import Column, Schema
from pandas_schema.validation import *
from cpipe.design import Design
from pathlib import Path

class FileExistsValidation(CustomElementValidation):
    def __init__(self):
        super().__init__(lambda x: Path(x).exists(), 'is not an existing file')

class FieldDeprecatedValidation(CustomElementValidation):
    def __init__(self):
        super().__init__(lambda: False, 'this field is deprecated and all data contained in it will be ignored')

class EmailValidation(CustomElementValidation):
    def __init__(self):
        super().__init__(self.__class__.valid_email, 'is an invalid email address')

    @staticmethod
    def valid_email(email):
        parsed = parseaddr(email)
        return all(parsed) and '@' in parsed[1]

class FastqListValidation(CustomElementValidation):
    def __init__(self):
        super().__init__(self.is_fastq_list, 'is not a list of existing files with the .fastq.gz suffix')

    def is_fastq_list(self, list):
        split = list.split(',')
        for fastq in split:
            fastq_path = Path(fastq)
            if fastq_path.exists() and set(fastq_path.suffixes).issubset({'.fastq', '.gz'}):
                continue
            else:
                return False

        return True

class PrioritisedGenesValidation(MatchesPatternValidation):
    def __init__(self):
        gene = r'[A-Z0-9]'
        priority = r'[1-4]:(?:{gene})(?:,{gene})*'.format(gene=gene)
        super().__init__(r'{priority}(?: {priority})'.format(priority=priority))

class MgDateListValidation(CustomElementValidation):
    def __init__(self):
        super().__init__(self.__class__.valid_date_list,
                         'is not a comma separated list of dates in the format "YYYYMMDD"')

    def valid_date_list(dates):
        split = dates.split(',')
        for date in split:
            try:
                datetime.datetime.strptime(date, '%Y%m%d')
            except:
                return False

        return True

cohorts = [design.name for design in Design.list_all()]

def get_schema(is_mgha):
    schema = Schema([
        Column('Batch', [InverseValidation(MatchesPatternValidation('_'), message='underscores are not allowed in batch IDs!')]),
        Column('Sample_ID', [InverseValidation(MatchesPatternValidation('_'), message='underscores are not allowed in sample IDs!')]),
        Column('DNA_Tube_ID', allow_empty=True),
        Column('Sex', [InListValidation(['male', 'female', 'other'], case_sensitive=False)]),
        Column('DNA_Concentration', [CanConvertValidation(float)], allow_empty=True),
        Column('DNA_Volume', [CanConvertValidation(float)], allow_empty=True),
        Column('DNA_Quantity', [CanConvertValidation(float)], allow_empty=True),
        Column('DNA_Quality', [CanConvertValidation(float)], allow_empty=True),
        Column('DNA_Date', [MgDateListValidation()], allow_empty=True),
        Column('Cohort', [InListValidation(cohorts)]),
        Column('Sample_Type',
               [InListValidation(['Normal', 'Tumour'])]
               if is_mgha
               else [], allow_empty=True),
        Column('Fastq_Files', [FastqListValidation()]),
        Column('Prioritised_Genes', [PrioritisedGenesValidation()], allow_empty=True),
        Column('Consanguinity',
               [InListValidation(['no', 'yes', 'suspected', 'unknown'], case_sensitive=False)],
               allow_empty=True
               ),
        Column('Variants_File', [FileExistsValidation()], allow_empty=True),
        Column('Pedigree_File', [InListValidation(['import', 'exclude']) | MatchesPatternValidation(r'\w+=\w+,\w+')],
               allow_empty=True),
        Column('Ethnicity', [InListValidation(['Unknown', 'European', 'African', 'Asian'])], allow_empty=True),
        Column('VariantCall_Group', [MatchesPatternValidation('[\w\d](?:,[\w\d])*')], allow_empty=True),
        Column('Capture_Date', [MgDateListValidation()], allow_empty=True),
        Column('Sequencing_Date', [MgDateListValidation()], allow_empty=True),
        Column('Mean_Coverage', [CanConvertValidation(float)], allow_empty=True),
        Column('Duplicate_Percentage', [CanConvertValidation(float)], allow_empty=True),
        Column('Machine_ID', allow_empty=True),
        Column('DNA_Extraction_Lab',
               [InListValidation(['RMH', 'MH', 'PMCC', 'AH', 'VCGS', 'CTP', 'Coriell'])]
               if is_mgha
               else [], allow_empty=True),

        Column('Sequencing_Lab',
               [InListValidation(['PMCC', 'AGRF', 'MH', 'VCGS'])]
               if is_mgha
               else [], allow_empty=True),
        Column('Exome_Capture', allow_empty=True),
        Column('Library_Preparation', allow_empty=True),
        Column('Barcode_Pool_Size', [CanConvertValidation(int)], allow_empty=True),
        Column('Read_Type', [MatchesPatternValidation('\d+(?:SE|PE)', message='The read type must consist of the read '
                                                                            'length in bp, followed by "SE" or "PE" '
                                                                            'for single of paired end reads, e.g. '
                                                                            '100PE')],
               allow_empty=True),
        Column('Machine_Type', allow_empty=True),
        Column('Sequencing_Chemistry', allow_empty=True),
        Column('Sequencing_Software', allow_empty=True),
        Column('Demultiplex_Software', allow_empty=True),
        Column('Hospital_Centre',
               [InListValidation(['RMH', 'PMCC', 'RCH', 'AH', 'MH'])]
               if is_mgha
               else [], allow_empty=True),

        Column('Sequencing_Contact', [EmailValidation()], allow_empty=True),
        Column('Pipeline_Contact', [EmailValidation()], allow_empty=True),
        Column('Notes', allow_empty=True),
    ])

    for column in schema.columns:
        column.validations.extend([LeadingWhitespaceValidation(), TrailingWhitespaceValidation()])

    return schema
