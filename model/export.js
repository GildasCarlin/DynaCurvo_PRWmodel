/**
 * File              : export.js
 * Original article  : doi
 * Authors           : Gildas Carlin, Ian Manifacier and Jean-Louis Milan
 * Date              : 05/05/2025
 * Description       : Simulation data is exported in a .txt file downloaded in the Downloads folder.
 */
//

let exported = false;
function data_out(){
  if (!exported && t>iTmax){
    exported = true;
    saveStrings(table_data, file_name_data);
  }
}

function export_the_data(){
  if (data){
    data_out(table_data);
  }
}


