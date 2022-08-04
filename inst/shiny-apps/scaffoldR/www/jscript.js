

$(document).ready(function(){
  
// confirm befors closing the app
window.onbeforeunload = function(e) {
  return '';
}; 

// show/hide plots and table  
$('#dtInput').click(function(){
  $('#DataInput').slideToggle();
});

 
$('#dtBin').click(function(){
  $('#DtBin').slideToggle();
});

$('#faInput').click(function(){
  $('#FAInput').slideToggle();
});

$('#vwData').click(function(){
  $('#DataView').slideToggle();
});

$('#vwData1').click(function(){
  $('#DataView1').slideToggle();
});

$('#vwStats').click(function(){
  $('#VwLnkData').slideUp();
  $('#VwStats').slideToggle();
});

$('#vwLnkData').click(function(){
   $('#VwStats').slideUp();
  $('#VwLnkData').slideToggle();
  
});

$('#intConfig1').click(function(){
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
  $('#VwJointMap').slideToggle();
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
});

$('#scafConfig').click(function(){
  $('#ScafConfig').slideToggle();
});

  // reset the app when the button pressed
  $('#reset').click(function(){
  history.go(0);
});


// tooltips
$("#dir1").attr('title', 'Scaffolding direction');
$("#mapSrc").attr('title', 'Data source');
$("#exitZoom1").attr('title', 'Exit zoom');
$("#exitZoom2").attr('title', 'Exit zoom');
$("#clearBrush1").attr('title', 'Clear break line');
$("#clearBrush2").attr('title', 'Clear break line');
$("#cut1").attr('title', 'Break sequence');
$("#cut2").attr('title', 'Break sequence');
$("#svChanges").attr('title', 'Save change to the set data directory');
$("#export").attr('title', 'Save to data directory');
$("#edit").attr('title', 'Edit');
$("#clipbtn").attr('title', 'Copy to clipboard');
$("#erase").attr('title', 'Clear list');
$("#reverse").attr('title', 'Reverse scaffold');
$("#check").attr('title', 'Ok');
});

