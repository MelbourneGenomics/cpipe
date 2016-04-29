/////////////////////////////////////////////////////////////////////////////////
//
// This file is part of Cpipe.
//
// Cpipe is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, under version 3 of the License, subject
// to additional terms compatible with the GNU General Public License version 3,
// specified in the LICENSE file that is part of the Cpipe distribution.
//
// Cpipe is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with Cpipe.  If not, see <http://www.gnu.org/licenses/>.
//
/////////////////////////////////////////////////////////////////////////////////
//
// Cpipe Main Pipeline Script
//
/////////////////////////////////////////////////////////////////////////////////

// remove spaces from gene lists and point to a new sample metadata file
// note that this isn't run through bpipe
correct_sample_metadata_file = {
    def target = new File('results')
    if( !target.exists() ) {
        target.mkdirs()
    }
    [ "sh", "-c", "python $SCRIPTS/correct_sample_metadata_file.py < $it > results/samples.corrected" ].execute().waitFor()
    return "results/samples.corrected"
}

/////////////////////////////////////////////////////////
// helper functions
/////////////////////////////////////////////////////////

// classify samples as singleton, trio, or member of trio
List find_sample_types(sample_info) {
    println "finding samples"

    // all samples
    all_samples = sample_info.keySet()

    // proband samples
    trio_samples = []
    proband_samples = all_samples.findAll { 
        if (sample_info[it].pedigree != "") {
            new_members = sample_info[it].pedigree.tokenize(';')[0].tokenize('=')[1].tokenize(','); // fid=na12877,na12878
            trio_samples.addAll(new_members); // add to members
            return true;
        }
        else {
            return false;
        }
    }

    individual_samples = all_samples.collect()
    individual_samples.removeAll(trio_samples) // includes probands
    println "done finding samples"
    return [ all_samples, proband_samples, trio_samples, individual_samples ]
}

// log changes to stage status
void stage_status(stage_name, stage_status, sample) {
    String current = new Date().format("yyMMdd-HHmmss")
    println("${current}: ${stage_name}: ${stage_status} (${sample})")
}
