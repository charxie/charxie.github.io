(this.webpackJsonpaladdin=this.webpackJsonpaladdin||[]).push([[0],{42:function(e,t,n){},43:function(e,t,n){},50:function(e,t,n){"use strict";n.r(t);var r=n(8),a=n.n(r),c=n(34),i=n.n(c),o=(n(42),n(43),n(15)),s=n(37),l=n(24),j=n(30),u=n(29),h=n(1);Object(u.b)();var d=Object(l.a)(Object(j.devtools)(Object(j.persist)((function(e,t,n){var r=function(t){return e(Object(u.a)(t))};return{set:r,worlds:{},getWorld:function(e){return t().worlds[e]},createNewWorld:function(){r((function(e){var t=[];t.push({type:"Foundation",cx:0,cy:0,lx:2,ly:4,height:.1,id:"f1"}),t.push({type:"Foundation",cx:1,cy:2,lx:2,ly:2,height:.2,id:"f2"});var n={name:"default",elements:t,cameraPosition:new h.Vector3(0,0,5)};e.worlds[n.name]=n}))}}}),{name:"aladdin-storage"}))),b=n(5);Object(o.b)({OrbitControls:s.a});var O=function(){var e=d((function(e){return e.set})),t=Object(o.e)(),n=t.camera,a=t.gl.domElement,c=Object(r.useRef)(null);Object(r.useEffect)((function(){return c.current&&(c.current.target.set(0,0,0),c.current.addEventListener("end",i)),function(){var e;null===(e=c.current)||void 0===e||e.removeEventListener("end",i)}}),[]);var i=function(){e((function(e){var t=e.worlds.default;t&&(t.cameraPosition.x=n.position.x,t.cameraPosition.y=n.position.y,t.cameraPosition.z=n.position.z)}))};return Object(b.jsx)("orbitControls",{ref:c,args:[n,a],enableZoom:!0,maxAzimuthAngle:Math.PI,minAzimuthAngle:-Math.PI})},f=n(10),x=n(13),p=n.p+"static/media/daysky.66925efd.jpg",m=n.p+"static/media/nightsky.2ee37da6.jpg",g=function(e){var t=e.type,n=void 0===t?"day sky":t,a=Object(x.a)(e,["type"]),c=Object(r.useRef)(null),i=Object(r.useMemo)((function(){var e,t=new h.TextureLoader;switch(n){case"night sky":e=t.load(m);break;default:e=t.load(p)}return e}),[n]);return Object(b.jsxs)("mesh",Object(f.a)(Object(f.a)({},a),{},{ref:c,name:"Sky",scale:1,onClick:function(e){e.intersections.length>0&&(e.intersections[0].object===c.current&&console.log("Sky clicked"))},children:[Object(b.jsx)("sphereGeometry",{args:[1e3,16,16,0,2*Math.PI,0,Math.PI/2+.01]}),Object(b.jsx)("meshBasicMaterial",{map:i,side:h.DoubleSide,opacity:1,color:"skyblue"})]}))},w=function(e){var t=e.endPoint,n=void 0===t?1e3:t,a=Object(x.a)(e,["endPoint"]),c=Object(r.useRef)(null),i=[];i.push(new h.Vector3(-n,0,0)),i.push(new h.Vector3(n,0,0));var o=(new h.BufferGeometry).setFromPoints(i),s=new h.LineBasicMaterial({color:"red",linewidth:10}),l=[];l.push(new h.Vector3(0,-n,0)),l.push(new h.Vector3(0,n,0));var j=(new h.BufferGeometry).setFromPoints(l),u=new h.LineBasicMaterial({color:"green",linewidth:10}),d=[];d.push(new h.Vector3(0,0,-n)),d.push(new h.Vector3(0,0,n));var O=(new h.BufferGeometry).setFromPoints(d),p=new h.LineBasicMaterial({color:"blue",linewidth:10});return Object(b.jsxs)("mesh",Object(f.a)(Object(f.a)({},a),{},{ref:c,children:[Object(b.jsx)("lineSegments",{args:[o,s]}),Object(b.jsx)("lineSegments",{args:[j,u]}),Object(b.jsx)("lineSegments",{args:[O,p]})]}))},v=n(38),y=function(e){var t=e.scale,n=void 0===t?.01:t,a=Object(x.a)(e,["scale"]),c=Object(o.d)(v.a,"static/assets/compass.obj"),i=Object(o.d)(h.FontLoader,"static/fonts/helvetiker_regular.typeface.json"),s=Object(r.useRef)(null),l=Object(o.e)().camera;Object(o.c)((function(e){if(s.current){var t=new h.Vector3(.88,-.8,0).unproject(l);s.current.position.set(t.x,t.y,t.z)}}));var j={font:i,height:0,size:.005},u=new h.MeshBasicMaterial({color:"white"}),d=new h.MeshBasicMaterial({color:"red"});return Object(b.jsxs)("mesh",Object(f.a)(Object(f.a)({},a),{},{ref:s,rotation:new h.Euler(-Math.PI/2,0,0),children:[Object(b.jsx)("mesh",{position:[-.001,.02,0],material:u,children:Object(b.jsx)("textGeometry",{args:["N",j]})}),Object(b.jsx)("mesh",{position:[-.0015,-.025,0],material:u,children:Object(b.jsx)("textGeometry",{args:["S",j]})}),Object(b.jsx)("mesh",{position:[-.025,-.002,0],material:u,children:Object(b.jsx)("textGeometry",{args:["W",j]})}),Object(b.jsx)("mesh",{position:[.02,-.002,0],material:u,children:Object(b.jsx)("textGeometry",{args:["E",j]})}),Object(b.jsx)("primitive",{object:c,scale:n,material:d})]}))},M=n(14),P=n(52),z=n(53),S=function(e){var t=e.id,n=e.cx,a=e.cy,c=e.lx,i=void 0===c?1:c,o=e.ly,s=void 0===o?1:o,l=e.height,j=void 0===l?.1:l,u=e.color,O=void 0===u?"gray":u,f=e.lineColor,x=void 0===f?"black":f,p=e.hovered,m=void 0!==p&&p,g=e.selected,w=void 0!==g&&g,v=d((function(e){return e.set})),y=Object(r.useRef)(),S=Object(r.useRef)(),I=Object(r.useRef)(),T=Object(r.useRef)(),V=Object(r.useRef)(),C=new h.Vector3(n-i/2,j/2,a-s/2),_=new h.Vector3(n-i/2,j/2,a+s/2),k=new h.Vector3(n+i/2,j/2,a-s/2),B=new h.Vector3(n+i/2,j/2,a+s/2),E=.002,R=function(e){v((function(n){var r=n.worlds.default;if(r){var a,c=Object(M.a)(r.elements);try{for(c.s();!(a=c.n()).done;){var i=a.value;if(i.id===t){i.hovered=e;break}}}catch(o){c.e(o)}finally{c.f()}}}))};return Object(b.jsxs)("group",{children:[Object(b.jsx)(P.a,{castShadow:!0,receiveShadow:!0,ref:y,name:"Foundation",onClick:function(e){e.intersections.length>0&&(e.intersections[0].object===y.current&&v((function(e){var n=e.worlds.default;if(n){var r,a=Object(M.a)(n.elements);try{for(a.s();!(r=a.n()).done;){var c=r.value;c.selected=c.id===t}}catch(i){a.e(i)}finally{a.f()}}})))},onPointerOver:function(e){e.intersections.length>0&&(e.intersections[0].object===y.current&&R(!0))},onPointerOut:function(e){R(!1)},args:[i,j,s],position:[n,j/2,a],children:Object(b.jsx)("meshStandardMaterial",{attach:"material",color:m?"lightGray":O})}),Object(b.jsxs)(b.Fragment,{children:[Object(b.jsx)(z.a,{points:[[C.x,j,C.z],[k.x,j,k.z]],children:Object(b.jsx)("lineBasicMaterial",{color:x})}),Object(b.jsx)(z.a,{points:[[k.x,j,k.z],[B.x,j,B.z]],children:Object(b.jsx)("lineBasicMaterial",{color:x})}),Object(b.jsx)(z.a,{points:[[B.x,j,B.z],[_.x,j,_.z]],children:Object(b.jsx)("lineBasicMaterial",{color:x})}),Object(b.jsx)(z.a,{points:[[_.x,j,_.z],[C.x,j,C.z]],children:Object(b.jsx)("lineBasicMaterial",{color:x})}),Object(b.jsx)(z.a,{points:[[C.x,E,C.z],[k.x,E,k.z]],children:Object(b.jsx)("lineBasicMaterial",{color:x})}),Object(b.jsx)(z.a,{points:[[k.x,E,k.z],[B.x,E,B.z]],children:Object(b.jsx)("lineBasicMaterial",{color:x})}),Object(b.jsx)(z.a,{points:[[B.x,E,B.z],[_.x,E,_.z]],children:Object(b.jsx)("lineBasicMaterial",{color:x})}),Object(b.jsx)(z.a,{points:[[_.x,E,_.z],[C.x,E,C.z]],children:Object(b.jsx)("lineBasicMaterial",{color:x})}),Object(b.jsx)(z.a,{points:[[C.x,E,C.z],[C.x,j,C.z]],children:Object(b.jsx)("lineBasicMaterial",{color:x})}),Object(b.jsx)(z.a,{points:[[k.x,E,k.z],[k.x,j,k.z]],children:Object(b.jsx)("lineBasicMaterial",{color:x})}),Object(b.jsx)(z.a,{points:[[_.x,E,_.z],[_.x,j,_.z]],children:Object(b.jsx)("lineBasicMaterial",{color:x})}),Object(b.jsx)(z.a,{points:[[B.x,E,B.z],[B.x,j,B.z]],children:Object(b.jsx)("lineBasicMaterial",{color:x})})]}),w&&Object(b.jsxs)(b.Fragment,{children:[Object(b.jsx)(P.c,{ref:S,args:[.1,6,6],position:C,children:Object(b.jsx)("meshStandardMaterial",{attach:"material",color:"white"})}),Object(b.jsx)(P.c,{ref:I,args:[.1,6,6],position:_,children:Object(b.jsx)("meshStandardMaterial",{attach:"material",color:"white"})}),Object(b.jsx)(P.c,{ref:T,args:[.1,6,6],position:k,children:Object(b.jsx)("meshStandardMaterial",{attach:"material",color:"white"})}),Object(b.jsx)(P.c,{ref:V,args:[.1,6,6],position:B,children:Object(b.jsx)("meshStandardMaterial",{attach:"material",color:"white"})})]})]})},I=function(e){var t=e.world;return Object(b.jsxs)("group",{children:[Object(b.jsx)(P.a,{castShadow:!0,receiveShadow:!0,args:[1,2,1],position:[0,1,0],children:Object(b.jsx)("meshStandardMaterial",{attach:"material",color:"white"})}),t.elements.filter((function(e){return"Foundation"===e.type})).map((function(e){return Object(b.jsx)(S,Object(f.a)({},e),e.id)}))]})},T=function(e){var t=e.color,n=void 0===t?"forestgreen":t,a=(Object(x.a)(e,["color"]),d((function(e){return e.set}))),c=Object(r.useRef)();return Object(b.jsx)(P.b,{receiveShadow:!0,ref:c,name:"Ground",onClick:function(e){e.intersections.length>0&&(e.intersections[0].object===c.current&&a((function(e){var t=e.worlds.default;if(t){var n,r=Object(M.a)(t.elements);try{for(r.s();!(n=r.n()).done;)n.value.selected=!1}catch(a){r.e(a)}finally{r.f()}}})))},rotation:[-Math.PI/2,0,0],position:[0,0,1],args:[1e4,1e4],children:Object(b.jsx)("meshStandardMaterial",{side:h.DoubleSide,attach:"material",color:n})})},V=n(12),C=n(0),_=n(4),k=function(){function e(){Object(C.a)(this,e)}return Object(_.a)(e,null,[{key:"UNIT_VECTOR_POS_X",get:function(){return new h.Vector3(1,0,0)}},{key:"UNIT_VECTOR_NEG_X",get:function(){return new h.Vector3(-1,0,0)}},{key:"UNIT_VECTOR_POS_Y",get:function(){return new h.Vector3(0,1,0)}},{key:"UNIT_VECTOR_NEG_Y",get:function(){return new h.Vector3(0,-1,0)}},{key:"UNIT_VECTOR_POS_Z",get:function(){return new h.Vector3(0,0,1)}},{key:"UNIT_VECTOR_NEG_Z",get:function(){return new h.Vector3(0,0,-1)}},{key:"ZERO_TOLERANCE",get:function(){return 1e-4}},{key:"HALF_PI",get:function(){return Math.PI/2}},{key:"TWO_PI",get:function(){return 2*Math.PI}},{key:"toRadians",value:function(e){return e*(Math.PI/180)}},{key:"toDegrees",value:function(e){return e*(180/Math.PI)}},{key:"sphericalToCartesianZ",value:function(e){var t=e.x*Math.cos(e.z),n=t*Math.cos(e.y),r=t*Math.sin(e.y),a=e.x*Math.sin(e.z);return e.set(n,r,a),e}}]),e}(),B=23.45/180*Math.PI,E=function(e,t,n){var r=Math.asin(Math.sin(t)*Math.sin(n)+Math.cos(t)*Math.cos(e)*Math.cos(n)),a=Math.sin(e)*Math.cos(t),c=Math.cos(n)*Math.sin(t)-Math.cos(e)*Math.cos(t)*Math.sin(n),i=Math.atan2(c,a),o=new h.Vector3(5,i,r);k.sphericalToCartesianZ(o),o.setX(-o.x);var s=o.z;return o.z=o.y,o.y=s,o},R=function(e){var t=e.date,n=(void 0===t&&new Date,e.latitude),a=void 0===n?42/180*Math.PI:n,c=(Object(x.a)(e,["date","latitude"]),Object(r.useState)(0)),i=Object(V.a)(c,2),o=i[0],s=(i[1],Object(r.useState)(0)),l=Object(V.a)(s,2),j=l[0],u=(l[1],Object(r.useState)(new h.Vector3)),d=Object(V.a)(u,2),O=d[0],f=d[1];Object(r.useEffect)((function(){return f(E(j,o,a)),function(){}}),[]);var p=Object(r.useMemo)((function(){for(var e=[],t=[],n=2*Math.PI/72,r=0,a=0;a<k.TWO_PI+n/2;a+=n){var c=Math.min(a,k.TWO_PI),i=.3;e.push(k.sphericalToCartesianZ(new h.Vector3(5,c,0))),e.push(k.sphericalToCartesianZ(new h.Vector3(5+i,c,0))),e.push(k.sphericalToCartesianZ(new h.Vector3(5,c+n,0))),e.push(k.sphericalToCartesianZ(new h.Vector3(5+i,c,0))),e.push(k.sphericalToCartesianZ(new h.Vector3(5+i,c+n,0))),e.push(k.sphericalToCartesianZ(new h.Vector3(5,c+n,0)));var o=void 0;k.TWO_PI-c>k.ZERO_TOLERANCE&&(i=r%3===0?.5:.3,(o=new h.Vector3(5,c,0)).z=.002,t.push(k.sphericalToCartesianZ(o)),(o=new h.Vector3(5+i,c,0)).z=.002,t.push(k.sphericalToCartesianZ(o))),r++}for(var s=3*e.length,l=new Float32Array(s),j=new Float32Array(s),u=new Float32Array(s),d=0;d<e.length;d++){var b=3*d;l[b]=e[d].x,l[b+1]=e[d].y,l[b+2]=e[d].z,j[b]=0,j[b+1]=1,j[b+2]=0;var O=Math.floor(d/18)%2===0?.2:1;u[b]=O,u[b+1]=O,u[b+2]=O}return[l,j,u,t]}),[]),m=Object(V.a)(p,4),g=m[0],w=m[1],v=m[2],y=m[3],M=Object(r.useMemo)((function(){for(var e=k.TWO_PI/96,t=[],n=-Math.PI;n<Math.PI+e/2;n+=e){var r=E(n,o,a);r.z>-.3&&t.push(r)}return t}),[]),P=Object(r.useMemo)((function(){for(var e=2*B/12,t=k.TWO_PI/96,n=new h.BufferGeometry,r=0,c=[],i=[],o=-B;o<B-e/2;o+=e)for(var s=-Math.PI;s<Math.PI-t/2;s+=t){var l=s+t,j=o+e;l>Math.PI&&(l=Math.PI),j>B&&(j=B);var u=E(s,o,a),d=E(l,o,a),b=E(l,j,a),O=E(s,j,a);(u.z>=0||d.z>=0||b.z>=0||O.z>=0)&&(c.push(u,d,b,O),i.push(r),i.push(r+1),i.push(r+2),i.push(r),i.push(r+2),i.push(r+3),r+=4)}return n.setFromPoints(c),n.setIndex(new h.BufferAttribute(new Uint16Array(i),1)),n}),[]);return Object(b.jsxs)("mesh",{rotation:new h.Euler(-Math.PI/2,0,0),children:[Object(b.jsxs)("mesh",{children:[Object(b.jsxs)("bufferGeometry",{attach:"geometry",children:[Object(b.jsx)("bufferAttribute",{attachObject:["attributes","position"],count:g.length/3,array:g,itemSize:3}),Object(b.jsx)("bufferAttribute",{attachObject:["attributes","normal"],count:w.length/3,array:w,itemSize:3}),Object(b.jsx)("bufferAttribute",{attachObject:["attributes","color"],count:v.length/3,array:v,itemSize:3})]}),Object(b.jsx)("meshBasicMaterial",{side:h.DoubleSide,vertexColors:!0,polygonOffset:!0,polygonOffsetFactor:-.7,polygonOffsetUnits:-2})]}),Object(b.jsx)("lineSegments",{args:[(new h.BufferGeometry).setFromPoints(y),new h.MeshBasicMaterial({color:0})]}),Object(b.jsxs)("mesh",{children:[Object(b.jsx)("lineSegments",{args:[(new h.BufferGeometry).setFromPoints(M),new h.MeshBasicMaterial({color:new h.Color(1,1,0),clippingPlanes:[new h.Plane(k.UNIT_VECTOR_POS_Y,0)]})]}),Object(b.jsx)("mesh",{args:[P,new h.MeshBasicMaterial({side:h.DoubleSide,color:new h.Color(1,1,0),transparent:!0,opacity:.5,clippingPlanes:[new h.Plane(k.UNIT_VECTOR_POS_Y,0)]})]}),Object(b.jsx)("mesh",{position:O,args:[new h.SphereGeometry(.25,20,20),new h.MeshBasicMaterial({color:4294967040})]})]})]})},F=function(){var e=d((function(e){return e.worlds})),t=d((function(e){return e.getWorld})),n=d((function(e){return e.createNewWorld})),a=e.default;Object(r.useEffect)((function(){t("default")||n()}),[]);var c=new h.Vector3(0,0,5);return a&&c.set(a.cameraPosition.x,a.cameraPosition.y,a.cameraPosition.z),console.log("x"),Object(b.jsxs)("div",{className:"App",children:[Object(b.jsxs)("div",{style:{backgroundColor:"lightblue",height:"60px",paddingTop:"10px",fontSize:"30px"},children:[Object(b.jsx)("img",{alt:"Logo",src:"static/assets/aladdin-logo.png",height:"50px",style:{verticalAlign:"middle"}}),Object(b.jsx)("span",{style:{paddingLeft:"20px",verticalAlign:"middle"},children:"Aladdin"})]}),Object(b.jsx)(o.a,{shadows:!0,camera:{position:c,fov:90},style:{height:"calc(100vh - 70px)",backgroundColor:"black"},children:Object(b.jsxs)(r.Suspense,{fallback:null,children:[Object(b.jsx)(O,{}),Object(b.jsx)("ambientLight",{intensity:.25}),Object(b.jsx)("directionalLight",{color:"white",position:[2,2,0],intensity:.5,castShadow:!0,"shadow-mapSize-height":512,"shadow-mapSize-width":512}),Object(b.jsx)("gridHelper",{args:[500,100,"gray","gray"]}),Object(b.jsx)(y,{}),Object(b.jsx)(w,{}),Object(b.jsx)(T,{}),Object(b.jsx)(g,{}),Object(b.jsx)(R,{date:new Date,latitude:42}),a&&Object(b.jsx)(I,{world:a})]})})]})},A=function(e){e&&e instanceof Function&&n.e(3).then(n.bind(null,54)).then((function(t){var n=t.getCLS,r=t.getFID,a=t.getFCP,c=t.getLCP,i=t.getTTFB;n(e),r(e),a(e),c(e),i(e)}))};i.a.render(Object(b.jsx)(a.a.StrictMode,{children:Object(b.jsx)(F,{})}),document.getElementById("root")),A()}},[[50,1,2]]]);
//# sourceMappingURL=main.840dd9cb.chunk.js.map