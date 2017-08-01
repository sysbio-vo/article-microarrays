// (C) Wolfgang Huber 2010-2011

// Script parameters - these are set up by R in the function 'writeReport' when copying the 
//   template for this script from arrayQualityMetrics/inst/scripts into the report.

var highlightInitial = [ false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false ];
var arrayMetadata    = [ [ "1", "515G.CEL", "PDM3", "cells derived from a 94-year-old donor", "94", "year", "Parkinson's disease", "male", "515G", "Parkinson's disease", "male" ], [ "2", "459F.CEL", "PDF1", "cells derived from a 81-year-old donor", "81", "year", "Parkinson's disease", "female", "459F", "Parkinson's disease", "female" ], [ "3", "515E.CEL", "PDM1", "cells derived from a 68-year-old donor", "68", "year", "Parkinson's disease", "male", "515E", "Parkinson's disease", "male" ], [ "4", "515H.CEL", "CF1", "cells derived from a 61-year-old donor", "61", "year", "normal", "female", "515H", "normal", "female" ], [ "5", "515J.CEL", "PDF3", "cells derived from a 69-year-old donor", "69", "year", "Parkinson's disease", "female", "515J", "Parkinson's disease", "female" ], [ "6", "515L.CEL", "CM3", "cells derived from a 62-year-old donor", "62", "year", "normal", "male", "515L", "normal", "male" ], [ "7", "459D.CEL", "CF3", "cells derived from a 74-year-old donor", "74", "year", "normal", "female", "459D", "normal", "female" ], [ "8", "459C.CEL", "CM1", "cells derived from a 61-year-old donor", "61", "year", "normal", "male", "515H", "normal", "male" ], [ "9", "515K.CEL", "PDF2", "cells derived from a 84-year-old donor", "84", "year", "Parkinson's disease", "female", "515K", "Parkinson's disease", "female" ], [ "10", "515D.CEL", "PDM2", "cells derived from a 67-year-old donor", "67", "year", "Parkinson's disease", "male", "515D", "Parkinson's disease", "male" ], [ "11", "515F.CEL", "PDF4", "cells derived from a 66-year-old donor", "66", "year", "Parkinson's disease", "female", "515F", "Parkinson's disease", "female" ], [ "12", "515B.CEL", "PDM4", "cells derived from a 78-year-old donor", "78", "year", "Parkinson's disease", "male", "515B", "Parkinson's disease", "male" ], [ "13", "515I.CEL", "CF2", "cells derived from a donor of unknown age", "not available", "  ", "normal", "female", "515I", "normal", "female" ], [ "14", "515C.CEL", "CM4", "cells derived from a 89-year-old donor", "89", "year", "normal", "male", "515C", "normal", "male" ], [ "15", "459B.CEL", "CM2", "cells derived from a 63-year-old donor", "63", "year", "normal", "male", "459B", "normal", "male" ], [ "16", "459A.CEL", "CF4", "cells derived from a 69-year-old donor", "69", "year", "normal", "female", "459A", "normal", "female" ] ];
var svgObjectNames   = [ "pca", "dens" ];

var cssText = ["stroke-width:1; stroke-opacity:0.4",
               "stroke-width:3; stroke-opacity:1" ];

// Global variables - these are set up below by 'reportinit'
var tables;             // array of all the associated ('tooltips') tables on the page
var checkboxes;         // the checkboxes
var ssrules;


function reportinit() 
{
 
    var a, i, status;

    /*--------find checkboxes and set them to start values------*/
    checkboxes = document.getElementsByName("ReportObjectCheckBoxes");
    if(checkboxes.length != highlightInitial.length)
	throw new Error("checkboxes.length=" + checkboxes.length + "  !=  "
                        + " highlightInitial.length="+ highlightInitial.length);
    
    /*--------find associated tables and cache their locations------*/
    tables = new Array(svgObjectNames.length);
    for(i=0; i<tables.length; i++) 
    {
        tables[i] = safeGetElementById("Tab:"+svgObjectNames[i]);
    }

    /*------- style sheet rules ---------*/
    var ss = document.styleSheets[0];
    ssrules = ss.cssRules ? ss.cssRules : ss.rules; 

    /*------- checkboxes[a] is (expected to be) of class HTMLInputElement ---*/
    for(a=0; a<checkboxes.length; a++)
    {
	checkboxes[a].checked = highlightInitial[a];
        status = checkboxes[a].checked; 
        setReportObj(a+1, status, false);
    }

}


function safeGetElementById(id)
{
    res = document.getElementById(id);
    if(res == null)
        throw new Error("Id '"+ id + "' not found.");
    return(res)
}

/*------------------------------------------------------------
   Highlighting of Report Objects 
 ---------------------------------------------------------------*/
function setReportObj(reportObjId, status, doTable)
{
    var i, j, plotObjIds, selector;

    if(doTable) {
	for(i=0; i<svgObjectNames.length; i++) {
	    showTipTable(i, reportObjId);
	} 
    }

    /* This works in Chrome 10, ssrules will be null; we use getElementsByClassName and loop over them */
    if(ssrules == null) {
	elements = document.getElementsByClassName("aqm" + reportObjId); 
	for(i=0; i<elements.length; i++) {
	    elements[i].style.cssText = cssText[0+status];
	}
    } else {
    /* This works in Firefox 4 */
	var success = false;
	i = 0; 
	/* Some of this looping could already be cached in reportInit() */
	while( (!success) & (i < ssrules.length) ) {
	    selector = ssrules[i].selectorText;  // The selector 
            if (!selector) 
		continue; // Skip @import and other nonstyle rules
            if (selector == (".aqm" + reportObjId)) {
		success = true; 
		ssrules[i].style.cssText = cssText[0+status];
	    } else {
		i++;
	    }
	}
    }

}

/*------------------------------------------------------------
   Display of the Metadata Table
  ------------------------------------------------------------*/
function showTipTable(tableIndex, reportObjId)
{
    var rows = tables[tableIndex].rows;
    var a = reportObjId - 1;

    if(rows.length != arrayMetadata[a].length)
	throw new Error("rows.length=" + rows.length+"  !=  arrayMetadata[array].length=" + arrayMetadata[a].length);

    for(i=0; i<rows.length; i++) 
 	rows[i].cells[1].innerHTML = arrayMetadata[a][i];
}

function hideTipTable(tableIndex)
{
    var rows = tables[tableIndex].rows;

    for(i=0; i<rows.length; i++) 
 	rows[i].cells[1].innerHTML = "";
}


/*------------------------------------------------------------
  From module 'name' (e.g. 'density'), find numeric index in the 
  'svgObjectNames' array.
  ------------------------------------------------------------*/
function getIndexFromName(name) 
{
    var i;
    for(i=0; i<svgObjectNames.length; i++)
        if(svgObjectNames[i] == name)
	    return i;

    throw new Error("Did not find '" + name + "'.");
}


/*------------------------------------------------------------
  SVG plot object callbacks
  ------------------------------------------------------------*/
function plotObjRespond(what, reportObjId, name)
{

    var a, i, status;

    switch(what) {
    case "show":
	i = getIndexFromName(name);
	showTipTable(i, reportObjId);
	break;
    case "hide":
	i = getIndexFromName(name);
	hideTipTable(i);
	break;
    case "click":
        a = reportObjId - 1;
	status = !checkboxes[a].checked;
	checkboxes[a].checked = status;
	setReportObj(reportObjId, status, true);
	break;
    default:
	throw new Error("Invalid 'what': "+what)
    }
}

/*------------------------------------------------------------
  checkboxes 'onchange' event
------------------------------------------------------------*/
function checkboxEvent(reportObjId)
{
    var a = reportObjId - 1;
    var status = checkboxes[a].checked;
    setReportObj(reportObjId, status, true);
}


/*------------------------------------------------------------
  toggle visibility
------------------------------------------------------------*/
function toggle(id){
  var head = safeGetElementById(id + "-h");
  var body = safeGetElementById(id + "-b");
  var hdtxt = head.innerHTML;
  var dsp;
  switch(body.style.display){
    case 'none':
      dsp = 'block';
      hdtxt = '-' + hdtxt.substr(1);
      break;
    case 'block':
      dsp = 'none';
      hdtxt = '+' + hdtxt.substr(1);
      break;
  }  
  body.style.display = dsp;
  head.innerHTML = hdtxt;
}
