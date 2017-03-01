<?php
    //**************************************//
    //                                      //
    //            Overlap finder            //
    //     Matthias Stahl, TU Muenchen      //
    //             Version 0.2              //
    //           18th April 2016            //
    //                                      //
    //**************************************//
    //
    //**************************************//
    // Description
    //**************************************//
    // Checks lists of proteomics data for the presence of unique keys
    //**************************************//

    // General variables

    // Specify file names to check

    /*$pgas_file_names = array(
        "test1.txt",
        "test2.txt"
    );*/

    // Definition of access key column
    $access_key = "Gene names";
    // $access_key = "Majority protein IDs";

    // Read PGAs collection
    $pgas_collection = array();
    foreach ($pgas_file_names as $pgas_file_name) {
        // Read list after list
        $pgas = array();
        $handle = fopen($pgas_file_name, "r");
        if ($handle) {
            $row_counter = 0;
            while (($data2 = fgetcsv($handle, 0, "\t")) !== FALSE) {
                $pgas[] = $data2;
                if ($row_counter == 0) {
                    // Save column headers
                    $column_headers[$pgas_file_name] = $pgas[0];
                    unset($pgas[0]);
                }
                $row_counter++;
            }
            fclose($handle);
        }
        // Add PGAs to collection
        $pgas_collection[$pgas_file_name] = $pgas;
        unset($pgas);
    }

    // Invert matrix
    $all_tables_inverted = array();
    foreach ($pgas_collection as $file_name => $pgas) {
        $invert_matrix = array();
        foreach ($pgas as $line_number => $line_content) {
            foreach ($line_content as $column_number => $column_content) {
                $invert_matrix[$line_number][$column_headers[$file_name][$column_number]] = $column_content;
            }
        }
        $all_tables_inverted[$file_name] = $invert_matrix;
    }

    // Comparison algorithm
    // Check all files provided
    foreach ($all_tables_inverted as $file_name => $lines) {
        $no_overlap = TRUE;
        $total_number_of_lines = count($lines);
        // Go through each line
        for ($basic_line_number = 1; $basic_line_number <= $total_number_of_lines; $basic_line_number++) {
            if (is_array($lines[$basic_line_number])) {
                $basic_entry = $lines[$basic_line_number];
                // Check specific line against all other lines
                for ($line_number = 1; $line_number <= $total_number_of_lines; $line_number++) {
                    if (is_array($lines[$line_number])) {
                        $entry = $lines[$line_number];
                        if ($basic_line_number != $line_number) {
                            if (($basic_entry[$access_key] != "") and ($entry[$access_key] != "")) { 
                                if (count(array_intersect(explode(";", $basic_entry[$access_key]), explode(";", $entry[$access_key]))) > 0) {
                                    $no_overlap = FALSE;
                                    echo("<strong>Overlap detected in</strong> ". $file_name ."<br />");
                                    echo("Key 1: ". $basic_entry[$access_key] ."<br />");
                                    echo("Key 2: ". $entry[$access_key] ."<br />");
                                    echo("Protein 1: ". $basic_entry["Protein names"] ."<br />");
                                    echo("Protein 2: ". $entry["Protein names"] ."<br />");
                                    echo("Line 1: ". $basic_line_number ."<br />");
                                    echo("Line 2: ". $line_number ."<br /><br />");
                                }
                            }
                        }
                    } else {
                        echo("No array!");
                    }
                }
            } else {
                echo("No array!");
            }
        }
        if ($no_overlap) {
            echo("<strong>There was no overlap found in</strong> ". $file_name ."<br />");
        }
        echo("<strong>Total number of lines was ". $total_number_of_lines ."</strong><br /><br />");
    }
?>