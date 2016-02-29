
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

germline_analysis_phase_1 = segment {
    // all samples do this
    variant_discovery

    // do this if an individual not in trio
    if (sample.is_not_in_trio()) {
        joint_call + 
        variant_annotation
    }
}

germline_analysis_phase_2 = segment {
    // does nothing
}
