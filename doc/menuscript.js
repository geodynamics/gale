
	var xV1 = 0;
	
	function stayHome(){	

		var nV = 0;
		if (!document.body.scrollTop){nV = document.documentElement.scrollTop}
		else {nV = document.body.scrollTop}
		document.getElementById('sidebar').style.top = xV1+nV+"px";
		setTimeout("stayHome()",50);
	}

	function init(){

		xV1 = document.getElementById('sidebar').offsetTop;
		stayHome();
	}

	window.onload=init;
	
