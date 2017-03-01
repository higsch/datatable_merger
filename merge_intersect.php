<?php
    //**************************************//
    //                                      //
    //           Merge_intersect            //
    //     Matthias Stahl, TU Muenchen      //
    //             Version 0.2              //
    //           4th April 2016             //
    //                                      //
    //**************************************//
    //
    //**************************************//
    // Description
    //**************************************//
    // Merges lists of proteomics data with unique keys
    //**************************************//

    // General variables

    // Specify file names here

    $pgas_file_names = array(
        "test1.txt",
        "test2.txt"
    );

    // Definition of access key column
    $access_key = "Majority protein IDs";
    //$access_key = "Fasta headers";
    //$access_key = "Gene names";
    //$access_key = "First ID wo Iso";
    //$access_key = "Fused protein IDs";

    // Output on screen and in csv file wanted?
    $write_output = TRUE;

    // If true, specify output file name
    $merge_file_name = "merge.csv";

    // Generate new inital key for new entries
    function generate_new_key($data) {
        $count = 0;
        foreach ($data as $matrix_lines) {
            if ($count < count($matrix_lines)) {
                $count = count($matrix_lines);
            }
        }
        $count += 1;
        return $count;
    }

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
    $counter = 0;
    $basic_matrix = array();
    $new_key = generate_new_key($all_tables_inverted);
    // Merge every matrices to the others still merged
    foreach ($all_tables_inverted as $file_name => $lines) {
        if ($counter == 0) {
            // First matrix
            $basic_matrix[$file_name] = $lines;
        } else {
            // 2nd+ matrices
            $tmp_matrix = array();
            // for each line in secondary matrix
            foreach ($lines as $line_number => $entry) {
                // Set exists flag
                $entry_exists = FALSE;
                // Generate array of search keys
                $secondary_entry_searcharray = explode(";", $entry[$access_key]);
                // for each column in basic matrix
                foreach ($basic_matrix as $basic_file_name => $basic_lines) {
                    // for each line in basic matrix
                    foreach ($basic_lines as $basic_line_number => $basic_entry) {
                        // Generate array of search keys
                        $basic_entry_searcharray = explode(";", $basic_entry[$access_key]);
                        // Main comparison of access keys
                        if (count(array_intersect($secondary_entry_searcharray, $basic_entry_searcharray)) > 0) {
                            // entry already exists
                            if (count(array_diff($tmp_matrix[$file_name][$basic_line_number], $entry)) > 0) {
                                echo("<strong>Entry already in array</strong><br />");
                                echo("Replicate match from ". $file_name ." to ". $basic_file_name ."<br />");
                                echo("Line number: ". $line_number ."/". $basic_line_number ."<br />");
                                echo("Gene name: ". $entry["Gene names"] ."<br />");
                                echo("Gene name (". $basic_file_name ."): ". $basic_entry["Gene names"] ."<br />");
                                echo("Former gene name: ". $tmp_matrix[$file_name][$basic_line_number]["Gene names"] ."<br /><br />");
                            }
                            // Assign entry to temporary matrix
                            $tmp_matrix[$file_name][$basic_line_number] = $entry;
                            $entry_exists = TRUE;
                        }
                    }
                }
                if (!$entry_exists) {
                    $tmp_matrix[$file_name][$new_key] = $entry;
                    $new_key++;
                }
            }
            // Assign temporary matrix to basic matrix
            foreach ($tmp_matrix as $file_name => $lines ){
                foreach ($lines as $line_number => $entry) {
                    $basic_matrix[$file_name][$line_number] = $entry;
                }
            }
            unset($tmp_matrix);
        }
        $counter++;
    }

    // Output in csv file
    // Create merge file
    if ($write_output) {
        $merge_csv = fopen($merge_file_name, "w");

        // Write column headers
        $line = array();
        foreach ($basic_matrix as $file_name => $content) {
            foreach ($column_headers[$file_name] as $column_name) {
                $line[] = $file_name ."_". $column_name;
            }
        }
        fputcsv($merge_csv, $line);

        // Write data
        $max_line = $new_key - 1;
        for ($i = 1; $i <= $max_line; $i++) {
            $line = array();
            foreach ($basic_matrix as $basic_file_name => $basic_lines) {
                if (array_key_exists($i, $basic_lines)) {
                    foreach ($column_headers[$basic_file_name] as $column_header) {
                        $line[] = $basic_lines[$i][$column_header];
                    }
                } else {
                    foreach ($column_headers[$basic_file_name] as $column_header) {
                        $line[] = "-";
                    }
                }
            }
            fputcsv($merge_csv, $line);
        }

        // Close file
        fclose($merge_csv);

        // Screen read-out
        echo("<h1>Merge_intersect v0.2</h1>");
        echo("<h3>developed in the Sieber research group, TUM, 2016</h3>");
        echo("<br />");
        echo("Access key for comparison: '". $access_key ."'<br />");
        echo("<pre>");
        print_r($basic_matrix);
        echo("</pre>");
    }
    echo("---FINISHED---")
?>