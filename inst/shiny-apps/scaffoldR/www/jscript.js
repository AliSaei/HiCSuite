

$(document).ready(function(){
  
// confirm befors closing the app
window.onbeforeunload = function(e) {
  return '';
};  
  
// show/hide plots and table  
$('#dtBin').click(function(){
  $('#DtBin').slideToggle();
});

$('#intConfig1').click(function(){
    $('#DtBin').slideUp();
  $('#IntConfig1').slideToggle();
});

$('#vwIntraMap').click(function(){
  $('#IntraMap').slideToggle();
	//$('#dSwitch').toggleClass('hide');
  //$('#Map1').toggleClass('col-sm-6 col-sm-4');
  //$('#Map2').toggleClass('col-sm-6 col-sm-8');
});

$('#vwMap1').click(function(){
  $('#IntraMap').slideUp();
  $('#Map1').slideToggle();
});

$('#vwScaf').click(function(){
  $('#IntraMap').slideUp();
  $('#VwScaf').slideToggle();
});


$('#vwJointMap').click(function(){
  $('#VwJointMap').toggle();
});

$('#vwIntData').click(function(){
  $('#VwIntData').slideToggle();
    $('#VwJointMap').slideUp();
});

$('#vwJointMap2').click(function(){
  $('#VwJointMap').slideUp();
   $('#VwIntData').slideUp();
  $('#VwJointMap2').slideToggle();
});


$('#intConfig2').click(function(){
  $('#IntConfig2').slideToggle();
  $('#ScafConfig').slideUp();
});

$('#scafConfig').click(function(){
  $('#IntConfig2').slideUp();
  $('#ScafConfig').slideToggle();
});

  // reset the app when the button pressed
  $('#reset').click(function(){
  history.go(0);
});


// tooltips
$("#dir1").attr('title', 'Scaffolding direction');
$("#mapSrc").attr('title', 'Data source');
$("#exitZoom").attr('title', 'Exit zoom');
$("#clearBrush").attr('title', 'Clear cut line');
$("#cut2").attr('title', 'Cut sequence');
$("#svChanges").attr('title', 'Save change to the set data directory');
$("#export").attr('title', 'Save to data directory');
$("#edit").attr('title', 'Edit');
$("#clipbtn").attr('title', 'Copy to clipboard');
$("#erase").attr('title', 'Clear list');
$("#check").attr('title', 'Ok');
});

