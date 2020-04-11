var DayCal = {
  DAYTIME: {'night': 'night', 'sunrise': 'sunrise', 'morning': 'morning', 'noon': 'noon', 'afternoon': 'afternoon', 'sunset': 'sunset', 'evening': 'evening', 'normal': 'normal'},
  theDate: new Date,
  theZone: 0,

  calDaytime: function(date, zone, glong, glat) {
    DayCal.theDate = date;
    DayCal.theZone = zone;
    mj = DayCal.mjd(date.getUTCDate(), date.getUTCMonth() + 1, date.getUTCFullYear(), 0.0);
    return DayCal.Cal(mj, zone, glong, glat);
  },

  mjd: function(day, month, year, hour) {
    var a, b;
    if (month <= 2) {
      month = month + 12;
      year = year - 1;
    }
    a = 10000.0 * year + 100.0 * month + day;
    if (a <= 15821004.1) {
      b = -2 * Math.floor((year + 4716)/4) - 1179;
    }
    else {
      b = Math.floor(year/400) - Math.floor(year/100) + Math.floor(year/4);
    }
    a = 365.0 * year - 679004.0;
    return (a + b + Math.floor(30.6001 * (month + 1)) + day + hour/24.0);
  },

  hrsmin: function(hours) {
    var hrs, h, m, dum;
    hrs = Math.floor(hours * 60 + 0.5)/ 60.0;
    h = Math.floor(hrs);
    m = Math.floor(60 * (hrs - h) + 0.5);
    if(h<10) h="0"+h;
    if(m<10) m="0"+m;
    dum = h +"时"+ m+"分";

    if (dum < 1000) dum = "0" + dum;
    if (dum <100) dum = "0" + dum;
    if (dum < 10) dum = "0" + dum;
    return dum;
  },

  ipart: function(x) {
    var a;
    if (x> 0) {
      a = Math.floor(x);
    }
    else {
      a = Math.ceil(x);
    }
    return a;
  },

  frac: function(x) {
    var a;
    a = x - Math.floor(x);
    if (a < 0) a += 1;
    return a;
  },

  range: function(x) {
    var a, b;
    b = x / 360;
    a = 360 * (b - DayCal.ipart(b));
    if (a  < 0 ) {
      a = a + 360
    }
    return a
  },

  quad: function(ym, yz, yp) {
    var nz, a, b, c, dis, dx, xe, ye, z1, z2, nz;
    var quadout = new Array;

    nz = 0;
    a = 0.5 * (ym + yp) - yz;
    b = 0.5 * (yp - ym);
    c = yz;
    xe = -b / (2 * a);
    ye = (a * xe + b) * xe + c;
    dis = b * b - 4.0 * a * c;
    if (dis > 0)	{
      dx = 0.5 * Math.sqrt(dis) / Math.abs(a);
      z1 = xe - dx;
      z2 = xe + dx;
      if (Math.abs(z1) <= 1.0) nz += 1;
      if (Math.abs(z2) <= 1.0) nz += 1;
      if (z1 < -1.0) z1 = z2;
    }
    quadout[0] = nz;
    quadout[1] = z1;
    quadout[2] = z2;
    quadout[3] = xe;
    quadout[4] = ye;
    return quadout;
  },

  lmst: function(mjd, glong) {
    var lst, t, d;
    d = mjd - 51544.5
    t = d / 36525.0;
    lst = DayCal.range(280.46061837 + 360.98564736629 * d + 0.000387933 *t*t - t*t*t / 38710000);
    return (lst/15.0 + glong/15);
  },

  minisun: function(t) {
    var p2 = 6.283185307, coseps = 0.91748, sineps = 0.39778;
    var L, M, DL, SL, X, Y, Z, RHO, ra, dec;
    var suneq = new Array;

    M = p2 * DayCal.frac(0.993133 + 99.997361 * t);
    DL = 6893.0 * Math.sin(M) + 72.0 * Math.sin(2 * M);
    L = p2 * DayCal.frac(0.7859453 + M / p2 + (6191.2 * t + DL)/1296000);
    SL = Math.sin(L);
    X = Math.cos(L);
    Y = coseps * SL;
    Z = sineps * SL;
    RHO = Math.sqrt(1 - Z * Z);
    dec = (360.0 / p2) * Math.atan(Z / RHO);
    ra = (48.0 / p2) * Math.atan(Y / (X + RHO));
    if (ra <0 ) ra += 24;
    suneq[1] = dec;
    suneq[2] = ra;
    return suneq;
  },

  sin_alt: function(iobj, mjd0, hour, glong, cglat, sglat) {
    var mjd, t, ra, dec, tau, salt, rads = 0.0174532925;
    var objpos = new Array;
    mjd = mjd0 + hour/24.0;
    t = (mjd - 51544.5) / 36525.0;
    if (iobj == 1) {
      objpos = minimoon(t);
    }
    else {
      objpos = DayCal.minisun(t);
    }
    ra = objpos[2];
    dec = objpos[1];
    tau = 15.0 * (DayCal.lmst(mjd, glong) - ra);
    salt = sglat * Math.sin(rads*dec) + cglat * Math.cos(rads*dec) * Math.cos(rads*tau);
    return salt;
  },

  getzttime: function(mjd, tz, glong) {
    var sglong, sglat, date, ym, yz, utrise, utset, j;
    var yp, nz, hour, z1, z2, iobj, rads = 0.0174532925;
    var quadout = new Array;

    sinho = Math.sin(rads * -0.833);
    date = mjd - tz/24;
    hour = 1.0;
    ym = DayCal.sin_alt(2, date, hour - 1.0, glong, 1, 0) - sinho;

    while(hour < 25) {
      yz = DayCal.sin_alt(2, date, hour, glong, 1, 0) - sinho;
      yp = DayCal.sin_alt(2, date, hour + 1.0, glong, 1, 0) - sinho;
      quadout = DayCal.quad(ym, yz, yp);
      nz = quadout[0];
      z1 = quadout[1];
      z2 = quadout[2];
      xe = quadout[3];
      ye = quadout[4];

      if (nz == 1) {
        if (ym < 0.0)
          utrise = hour + z1;
        else
          utset = hour + z1;

      }
      if (nz == 2) {
        if (ye < 0.0) {
          utrise = hour + z2;
          utset = hour + z1;
        }
        else {
          utrise = hour + z1;
          utset = hour + z2;
        }
      }
      ym = yp;
      hour += 2.0;
    }
    var zt=(utrise+utset)/2;
    if(zt<utrise)
      zt=(zt+12)%24;
    return zt;
  },

  Cal: function(mjd, tz, glong, glat) {
    var sglong, sglat, date, ym, yz, above, utrise, utset, j;
    var yp, nz, rise, sett, hour, z1, z2, iobj, rads = 0.0174532925;
    var quadout = [];

    sinho = Math.sin(rads * -0.833);
    sglat = Math.sin(rads * glat);
    cglat = Math.cos(rads * glat);
    date = mjd - tz/24;

    rise = false;
    sett = false;
    above = false;
    hour = 1.0;
    ym = DayCal.sin_alt(2, date, hour - 1.0, glong, cglat, sglat) - sinho;
    if (ym > 0.0) above = true;
    while(hour < 25 && (sett == false || rise == false)) {
      yz = DayCal.sin_alt(2, date, hour, glong, cglat, sglat) - sinho;
      yp = DayCal.sin_alt(2, date, hour + 1.0, glong, cglat, sglat) - sinho;
      quadout = DayCal.quad(ym, yz, yp);
      nz = quadout[0];
      z1 = quadout[1];
      z2 = quadout[2];
      xe = quadout[3];
      ye = quadout[4];

      if (nz == 1) {
        if (ym < 0.0) {
          utrise = hour + z1;
          rise = true;
        }
        else {
          utset = hour + z1;
          sett = true;
        }
      }

      if (nz == 2) {
        if (ye < 0.0) {
          utrise = hour + z2;
          utset = hour + z1;
        }
        else {
          utrise = hour + z1;
          utset = hour + z2;
        }
      }
      ym = yp;
      hour += 2.0;
    }
    if (rise == true || sett == true) {
      if(!rise || !sett) {
        return DayCal.DAYTIME.normal;
      }
      else {
        var h = DayCal.theDate.getUTCHours() + DayCal.theDate.getUTCMinutes() / 60;
        if(Math.abs(h - utrise) <= 1) {
          return DayCal.DAYTIME.sunrise;
        }
        else if(Math.abs(h - utset) <= 1) {
          return DayCal.DAYTIME.sunset;
        }
        else if(h > utrise && h < utset) {
          return DayCal.DAYTIME.noon;
        }
        else {
          return DayCal.DAYTIME.night;
        }
      }
    }
    else {
      if (above == true) {
        return DayCal.DAYTIME.noon;
      }
      else {
        return DayCal.DAYTIME.night;
      }
    }
  }
}

var Clock = {
  INITIATED: false,
  CANVAS_SIZE:[148, 148],
  UI:{},
  MYCLOCKS:[],
  BUTTON_COLORS:["btn-primary","btn-info","btn-success","btn-warning","btn-danger","btn-inverse"],
  UTC:new Date,
  DST_REGIONS:{},
  MONTHS:["January", "February", "March", "April", "May", "June", "July", "August", "September", "October", "November", "December"],
  HOURMODE:12,
  init:function () {
    if(Clock.INITIATED) {
      return;
    }
    Clock.INITIATED = true;
    this.ZONES = ZONES;
    this.DST_RULES = DST_RULES;
    this.CITIES = [];
    for (var e in this.ZONES) {
      this.CITIES.push(ZONES[e].name + " (" + ZONES[e].country + ")");
    }
    this.UTC = new Date;
    console.log('UTC:', this.UTC);
    this.initDST();
    this.loadMyClocks();
    this.initUI();
    this.initClocks();
    this.initDaytime();
  },
  initTime:function (e) {
  },
  initDST:function () {
    for (var e in this.ZONES) {
      if (this.ZONES[e].region) {
        var t = this.ZONES[e].region;
        var n = new Date(this.UTC.getTime() + this.ZONES[e].offset * 6e4);
        var r = this.DST_RULES[t]["dates"][n.getUTCFullYear()];
        var i = new Date(r[0] * 1e3), s = new Date(r[1] * 1e3);
        if (!this.DST_RULES[t].rev) {
          if (n > i && n < s) {
            this.ZONES[e].offset += 60;
          }
        } else {
          if (n > i && n < s) {
          } else {
            this.ZONES[e].offset += 60;
          }
        }
      }
    }
  },
  initUI:function () {
    this.UI.current_theme = "light";
    this.UI.themes = {dark:{circle:["#333333", "#111111"], ticks:"#fff", text:"#fff", arc:"#fff"}, light:{circle:["#333333", "#111111"], ticks:"#fff", text:"#fff", arc:"#000"}};
    if (navigator.userAgent.match(/phone/i) || navigator.userAgent.match(/android/i)) {
      $("body").addClass("phone");
    }
    this.UI.clocks = $("#clocks");
    this.UI.zone = $("#zone");
    this.UI.clockHash = {};
    this.UI.theme = function (e) {
      Clock.UI.current_theme = e;
      $("body").removeClass(e == "dark" ? "light" : "dark").addClass(e);
    };
    $("#add-form").submit(function () {
      var e = Clock.getIdFromName(Clock.UI.zone.val());
      Clock.UI.zone.val("");
      if (!e) {
        return false
      }
      if (Clock.MYCLOCKS.indexOf(e) != -1) {
        return false
      }
      Clock.addZone(e);
      Clock.createClock(e);
      Clock.renderClock(e);
      return false
    });
    $("#btn-reset").click(function () {
      localStorage.MYCLOCKS = "";
      localStorage.theme = "dark";
      document.location.href = "?reset=1"
    });
    $("#btn-theme").click(function () {
      var e = $(this).data("theme") == "dark" ? "light" : "dark";
      localStorage.theme = e;
      Clock.UI.theme(e);
      $(this).data("theme", e).html(e == "dark" ? "light" : "dark")
    });
    if (localStorage.theme) {
      this.UI.theme(localStorage.theme);
      $("#btn-theme").data("theme", localStorage.theme).html(localStorage.theme == "dark" ? "light" : "dark")
    } else {
      this.UI.theme("dark");
      $("#btn-theme").data("theme", "dark").html("light")
    }
    $("#btn-hourmode").click(function () {
      Clock.HOURMODE = $(this).data("mode") == 12 ? 24 : 12;
      $(this).data("mode", Clock.HOURMODE).html(Clock.HOURMODE + "h");
      localStorage.hourmode = Clock.HOURMODE
    });
    if (localStorage.hourmode) {
      this.HOURMODE = localStorage.hourmode;
      $("#btn-hourmode").data("mode", this.HOURMODE).html(this.HOURMODE + "h")
    } else {
      localStorage.hourmode = this.HOURMODE
    }
    this.UI.clocks.sortable({update:function () {
      Clock.saveReordered();
    }});
    this.UI.zone.focus(function () {
      if (this.value == "please enter a city") {
        this.value = ""
      }
    }).blur(function () {
        if (this.value == "") {
          this.value = "please enter a city"
        }
      })

    if(this.UI.zone.get(0)) {
      new IncrementalSearch(this.UI.zone.get(0), function (e, t) {
        if (t = t.toLowerCase()) {
          if (t.length > 2) {
            for (var n = -1, r = Clock.CITIES.length; ++n < r;) {
              for (var i = 0, s = []; i = Clock.CITIES[n].toLowerCase().indexOf(t, i) + 1; s[s.length] = i - 1);
              if (s.length)
                e.add(Clock.CITIES[n], s);
            }
          }
        }
        e.show();
      }, "autocomplete");
    }
  },
  initDaytime:function () {
    Clock.UTC = new Date(Clock.UTC.getTime() + 1e3);
    for (var t in Clock.MYCLOCKS) {
      Clock.renderDaytime(Clock.MYCLOCKS[t]);
    }

    /*window.setInterval(function () {
      Clock.UTC = new Date(Clock.UTC.getTime() + 1e3);
      for (var e in Clock.MYCLOCKS) {
        Clock.renderDaytime(Clock.MYCLOCKS[e]);
      }
    }, 18e2);*/
    Clock.UTC = new Date(Clock.UTC.getTime() + 1e3);
    for (var e in Clock.MYCLOCKS) {
      Clock.renderDaytime(Clock.MYCLOCKS[e]);
    }
  },
  initClocks:function () {
    Clock.UI.clocks.children().remove();
    for (var t in Clock.MYCLOCKS) {
      Clock.createClock(Clock.MYCLOCKS[t]);
      // render 1st time
      Clock.renderClock(Clock.MYCLOCKS[t]);
    }
    window.setInterval(function () {
      Clock.UTC = new Date(Clock.UTC.getTime() + 1e3);
      for (var e in Clock.MYCLOCKS) {
        Clock.renderClock(Clock.MYCLOCKS[e]);
      }
    }, 1e3);
  },
  removeClock:function (e) {
    Clock.removeZone(e);
    Clock.UI.clockHash[e].remove();
    delete Clock.UI.clockHash[e];
  },
  createClock:function (e) {
    var t = $('<Button class="btn btn-danger btn-small delete">').html('Delete').click(function () {
      Clock.removeClock(e);
      return false;
    });
    var n = $("<canvas>").width(this.CANVAS_SIZE[0]).height(this.CANVAS_SIZE[1]);
    var a = $('<div class="analog">').append(n);
    var r = $('<div class="clock">').attr('id', e).append($('<div class="LCD">').append($('<span class="label-city">').html(this.ZONES[e].name)).append($('<span class="label-utc">')).append($('<span class="label-time">')).append($('<span class="label-date">')).append($('<span class="label-weather">'))).data("id", e);
    r.append(t);
    this.UI.clocks.append(r.prepend(a));
    this.UI.clockHash[e] = r;
  },
  renderClock:function (e) {
    this.renderAnalog(e);
  },
  renderDaytime:function (e) {
    var t = new Date(this.UTC.getTime() + this.ZONES[e].offset * 6e4);
    var c = this.UI.clockHash[e];
    var r = DayCal.calDaytime(t, this.ZONES[e].offset / 60, this.ZONES[e].Za, this.ZONES[e].Ya);
    c.find('.analog').css('background', "url('/assets/images/face_" + r + ".png') no-repeat");
  },
  renderAnalog:function (e) {
    var t = new Date(this.UTC.getTime() + this.ZONES[e].offset * 6e4);
    var c = this.UI.clockHash[e];
    var n = this.UI.themes[this.UI.current_theme];
    var r = c.find('canvas').get(0).getContext("2d");
    var i = this.CANVAS_SIZE[0], s = this.CANVAS_SIZE[1], o = Math.min(i, s) / 2 - 10, u = i / 60;
    var l = this.ZONES[e].offset / 60;
    var a = t.getUTCHours(), f = "";
    if (this.HOURMODE == 12) {
      var f = a >= 12 ? " PM" : " AM", a = a > 12 ? a - 12 : a;
    }
    var h = t.getUTCMinutes(), p = t.getUTCSeconds();
    var d = t.getUTCHours();
    d = d >= 12 ? d - 12 : d;

    c.find('.label-utc').html("UTC " + (l > 0 ? "+" + l : l));
    c.find('.label-time').html(a.pad(2) + ":" + t.getUTCMinutes().pad(2) + f);
    c.find('.label-date').html(t.getUTCDate() + " " + this.MONTHS[t.getUTCMonth()].substring(0, 3));

    r.canvas.width = i;
    r.canvas.height = s;
    r.clearRect(0, 0, i, s);
    r.save();
    r.restore();
    r.translate(i / 2, s / 2);
    r.rotate(-Math.PI / 2);
    r.lineWidth = u;
    r.save();
    r.restore();
    r.strokeStyle = n.ticks;
    for (var c = 0; c < 12; c++) {
      r.beginPath();
      r.rotate(Math.PI / 6);
      r.moveTo(o - Math.ceil(o * .1), 0);
      r.lineTo(o, 0);
      r.stroke()
    }
    r.restore();
    r.save();
    r.rotate(Math.PI / 6 * d + Math.PI / 360 * h + Math.PI / 21600 * p);
    r.lineWidth = u;
    r.beginPath();
    r.moveTo(0, 0);
    r.lineTo(Math.ceil(o * .55), 0);
    r.stroke();
    r.restore();
    r.save();
    r.rotate(Math.PI / 30 * h + Math.PI / 1800 * p);
    r.beginPath();
    r.moveTo(0, 0);
    r.lineTo(Math.ceil(o * .8), 0);
    r.stroke();
    r.restore();
    r.save();
    r.rotate(p * Math.PI / 30);
    r.strokeStyle = "#ff0000";
    r.fillStyle = "#ff0000";
    r.lineWidth = 2;
    r.beginPath();
    r.moveTo(-Math.ceil(o * .2), 0);
    r.lineTo(o - Math.ceil(o * .1), 0);
    r.stroke();
    r.beginPath();
    r.arc(0, 0, Math.ceil(o * .06), 0, Math.PI * 2, true);
    r.fill();
    r.restore();
    r.save();
  },
  loadMyClocks:function () {
    this.MYCLOCKS = [];
    var e = [];
    if (localStorage.MYCLOCKS && JSON.parse(localStorage.MYCLOCKS).length > 0) {
      e = JSON.parse(localStorage.MYCLOCKS)
    }
    for (var t in e) {
      if (this.ZONES[e[t]]) {
        this.MYCLOCKS.push(e[t])
      }
    }
    this.commit();
    if (this.MYCLOCKS.length < 1) {
      $.each(["new_york.usa", "london.united_kingdom", "cairo.egypt", "new_delhi.india", "beijing.china", "tokyo.japan"], function (e, t) {
        Clock.addZone(t);
      });
    }
  },
  getIdFromName:function (e) {
    for (var t in ZONES) {
      if (ZONES[t].name + " (" + ZONES[t].country + ")" == e) {
        return t;
      }
    }
    return null;
  },
  addZone:function (e) {
    this.MYCLOCKS.push(e);
    this.commit();
  },
  removeZone:function (e) {
    this.MYCLOCKS.splice(this.MYCLOCKS.indexOf(e), 1);
    this.commit();
  },
  saveReordered:function () {
    this.MYCLOCKS = new Array;
    this.UI.clocks.find(".clock").each(function () {
      Clock.MYCLOCKS.push($(this).data("id"));
    });
    this.commit();
  },
  commit:function () {
    localStorage.MYCLOCKS = JSON.stringify(this.MYCLOCKS);
  }
};

Number.prototype.pad = function (e) {
  var t = "0" + this;
  return t.substr(t.length - e);
};

IncrementalSearch = function (e, t, n) {
  var r, i = this;
  (i.input = e).autocomplete = "off", i.callback = t || function () {
  }, i.className = n || "", i.hide(), i.visible = 0;
  for (r in {keydown:0, focus:0, blur:0, keyup:0, keypress:0}) {
    addEvent(e, r, i._handler, i);
  }
};

with ({p:IncrementalSearch.prototype}) {
  p.show = function () {
    for (var e = this, t = document.body.appendChild(e.c).style, n = e.input, r = n.offsetLeft, i = n.offsetTop + n.offsetHeight; n = n.offsetParent; r += n.offsetLeft, i += n.offsetTop);
    t.left = r + "px", t.top = i + "px", t.minWidth = "238px";
    e.l.length ? (t.display = "block", !e.visible && (e._callEvent("onshow"), ++e.visible), e.highlite(0)) : t.display = "none"
  };
  p.hide = function () {
    var e = this, t = (e.c && e.c.parentNode && e.c.parentNode.removeChild(e.c), e.c = document.createElement("div")).style;
    e.l = [], e.i = -1, e.c.className = e.className, t.position = "absolute", t.display = "none";
    e._old = null, e.visible && (e._callEvent("onhide"), --e.visible)
  };
  p.add = function (e, t, n) {
    var r = this, i = 0, s = document, o = r.l.length, u = r.input.value.length, a = (r.l[o] = [e, n, r.c.appendChild(s.createElement("div"))])[2];
    if (t instanceof Array || (t = [t]), a.i = o, a.className = "normal", !isNaN(t[0]))for (var f = -1, l = t.length; ++f < l; a.appendChild(s.createTextNode(e.substring(i, t[f]))).parentNode.appendChild(s.createElement("span")).appendChild(s.createTextNode(e.substring(t[f], i = t[f] + u))).parentNode.className = "highlited");
    for (t in a.appendChild(s.createTextNode(e.substr(i))), {click:0, mouseover:0})addEvent(a, t, r._handler, r)
  };
  p.highlite = function (e) {
    var t = this;
    t._invalid(e) || (t._invalid(t.i) || (t.l[t.i][2].className = "normal"), t.l[t.i = e][2].className += " selected", t._callEvent("onhighlite", t.l[e][0], t.l[e][1]))
  };
  p.select = function (e) {
    var t = this;
    t._invalid(e = isNaN(e) ? t.i : e) || (t._callEvent("onselect", t.input.value = t.l[t.i][0], t.l[e][1]), t.hide())
  };
  p.next = function () {
    var e = (e = this, e.highlite((e.i + 1) % e.l.length))
  };
  p.previous = function () {
    var e = (e = this, e.highlite((!e.i ? e.l.length : e.i) - 1))
  };
  p._fadeOut = function () {
    var e = (e = function () {
      arguments.callee.x.hide()
    }, e.x = this, setTimeout(e, 200))
  };
  p._handler = function (e) {
    var t = this, n = e.type, r = e.key;
    n == "focus" || n == "keyup" ? r != 40 && r != 38 && r != 13 && t._old != t.input.value && (t.hide(), t.callback(t, t.input.value)) : n == "keydown" ? r == 40 ? t.next() : r == 38 ? t.previous() : t._old = t.input.value : n == "keypress" ? r == 13 && t.select() : n == "blur" ? t._fadeOut() : n == "click" ? t.select() : t.highlite((/span/i.test((e = e.target).tagName) ? e.parentNode : e).i)
  };
  p._invalid = function (e) {
    return isNaN(e) || e < 0 || e >= this.l.length
  };
  p._callEvent = function (e) {
    var t = this;
    return t[e]instanceof Function ? t[e].apply(t, [].slice.call(arguments, 1)) : undefined;
  }
}

function addEvent(e, t, n, r) {
  var i = e[i = "_" + (t = "on" + t)] = e[i] || (e[t] ? [
    [e[t], e]
  ] : []), s, o, u;
  i[i.length] = [n, r || e], e[t] = function (t) {
    try {
      (t = t || event).preventDefault || (t.preventDefault = function () {
        t.returnValue = false
      });
      t.stopPropagation || (t.stopPropagation = function () {
        t.cancelBubble = true
      });
      t.target || (t.target = t.srcElement || null);
      t.key = (t.which + 1 || t.keyCode + 1) - 1 || 0;
    } catch (n) {
    }
    for (u = 1, n = i.length; n; i[--n] && (s = i[n][0], e = i[n][1], s.call ? o = s.call(e, t) : (e._ = s, o = e._(t), e._ = null), u &= o !== false));
    return t = null, !!u;
  }
}

Clock.init();
