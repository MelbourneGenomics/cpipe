#!/usr/bin/env groovy
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

CliBuilder cli = new CliBuilder()
cli.with {
    batch "batch to which samples belong", args:1, required: true
    disease "disease cohort to which samples belong", args:1, required: true
}

opts = cli.parse(args)
if(!opts)
        System.exit(0)

samples = SampleInfo.fromFiles(opts.arguments() as List)
samples.each { it.value.batch = opts.batch; it.value.target = opts.disease }

println(
    samples*.value*.toTsv().join("\n")
)
