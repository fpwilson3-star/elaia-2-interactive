const state = {
  table1: null,
  fig1: null,
  fig2: null,
  fig3: null,
  table2: null,
  table3: null,
};

const ARM_COLORS = {
  alert: "#1a5490",   // trial's "alert" arm
  usual: "#c1440e",   // trial's "usual care" arm
};

// ---------- utilities ----------

function escapeHtml(s) {
  return String(s)
    .replaceAll("&", "&amp;")
    .replaceAll("<", "&lt;")
    .replaceAll(">", "&gt;")
    .replaceAll('"', "&quot;");
}

function fmtPct(p) {
  return Number.isFinite(p) ? `${(100 * p).toFixed(1)}%` : "—";
}

function fmtP(p) {
  if (p == null) return "—";
  if (p < 0.0001) return p.toExponential(2);
  if (p < 0.01) return p.toFixed(3);
  return p.toFixed(2);
}

// Populate a <select> with <optgroup>s keyed by each entry's `group` field.
// `subgroups` is an object { key: {label, group, ...} } preserving insertion order.
function populateGroupedSelect(select, subgroups) {
  while (select.firstChild) select.removeChild(select.firstChild);
  const byGroup = new Map();
  for (const [key, sub] of Object.entries(subgroups)) {
    const g = sub.group || "Other";
    if (!byGroup.has(g)) byGroup.set(g, []);
    byGroup.get(g).push([key, sub]);
  }
  for (const [groupName, items] of byGroup) {
    const og = document.createElement("optgroup");
    og.label = groupName;
    for (const [key, sub] of items) {
      const opt = document.createElement("option");
      opt.value = key;
      opt.textContent = sub.label;
      og.appendChild(opt);
    }
    select.appendChild(og);
  }
}

async function fetchJSON(path) {
  const resp = await fetch(path);
  if (!resp.ok) throw new Error(`failed to load ${path}: ${resp.status}`);
  return resp.json();
}

// ---------- Table 1 ----------

function renderTable1(key) {
  const sub = state.table1.subgroups[key];
  document.getElementById("table1-n").textContent =
    `Alert: n = ${sub.n_alert.toLocaleString()} · Usual care: n = ${sub.n_usual.toLocaleString()}`;

  const table = document.createElement("table");
  table.className = "data-table";
  table.innerHTML = `
    <thead>
      <tr>
        <th>Characteristic</th>
        <th class="col-num">Alert (N = ${sub.n_alert.toLocaleString()})</th>
        <th class="col-num">Usual care (N = ${sub.n_usual.toLocaleString()})</th>
      </tr>
    </thead>`;
  const tbody = document.createElement("tbody");
  for (const section of sub.sections) {
    tbody.insertAdjacentHTML(
      "beforeend",
      `<tr class="section-row"><td colspan="3">${escapeHtml(section.title)}</td></tr>`
    );
    for (const row of section.rows) {
      tbody.insertAdjacentHTML(
        "beforeend",
        `<tr>
          <td>${escapeHtml(row.label)}</td>
          <td class="col-num">${escapeHtml(row.alert)}</td>
          <td class="col-num">${escapeHtml(row.usual)}</td>
        </tr>`
      );
    }
  }
  table.appendChild(tbody);
  document.getElementById("table1-container").replaceChildren(table);
}

function initTable1() {
  const select = document.getElementById("table1-subgroup");
  populateGroupedSelect(select, state.table1.subgroups);
  select.addEventListener("change", () => renderTable1(select.value));
  document.getElementById("table1-footnote").textContent = state.table1.footnote;
  renderTable1("overall");
}

// ---------- Figure plotting (shared by Fig 1 and Fig 2) ----------

function renderFigure(containerId, readoutId, nSummaryId, figData, cohortKey) {
  const cohort = figData.subgroups[cohortKey];
  const panels = cohort.panels;

  if (nSummaryId) {
    document.getElementById(nSummaryId).textContent =
      `Alert: n = ${cohort.n_alert.toLocaleString()} · Usual care: n = ${cohort.n_usual.toLocaleString()}`;
  }

  const xs = [];
  const alertY = [];
  const alertErr = [];
  const usualY = [];
  const usualErr = [];
  const alertText = [];
  const usualText = [];

  for (const panel of panels) {
    xs.push(panel.label);
    const a = panel.data.alert;
    const u = panel.data.usual;
    alertY.push(a.prop);
    alertErr.push([a.prop - a.ci_lo, a.ci_hi - a.prop]);
    usualY.push(u.prop);
    usualErr.push([u.prop - u.ci_lo, u.ci_hi - u.prop]);
    alertText.push(`n=${a.n}<br>events=${a.events}<br>${fmtPct(a.prop)} (${fmtPct(a.ci_lo)}–${fmtPct(a.ci_hi)})`);
    usualText.push(`n=${u.n}<br>events=${u.events}<br>${fmtPct(u.prop)} (${fmtPct(u.ci_lo)}–${fmtPct(u.ci_hi)})`);
  }

  const traces = [
    {
      name: "Alert",
      x: xs,
      y: alertY,
      type: "scatter",
      mode: "markers",
      marker: { size: 12, color: ARM_COLORS.alert, symbol: "circle" },
      error_y: {
        type: "data", symmetric: false,
        array: alertErr.map((e) => e[1]), arrayminus: alertErr.map((e) => e[0]),
        color: ARM_COLORS.alert, thickness: 2, width: 8,
      },
      hovertemplate: "<b>Alert</b><br>%{text}<extra></extra>",
      text: alertText,
      offsetgroup: "alert",
    },
    {
      name: "Usual care",
      x: xs,
      y: usualY,
      type: "scatter",
      mode: "markers",
      marker: { size: 12, color: ARM_COLORS.usual, symbol: "square" },
      error_y: {
        type: "data", symmetric: false,
        array: usualErr.map((e) => e[1]), arrayminus: usualErr.map((e) => e[0]),
        color: ARM_COLORS.usual, thickness: 2, width: 8,
      },
      hovertemplate: "<b>Usual care</b><br>%{text}<extra></extra>",
      text: usualText,
      offsetgroup: "usual",
    },
  ];

  const layout = {
    margin: { l: 60, r: 20, t: 10, b: 70 },
    yaxis: {
      title: { text: figData.outcome_label, standoff: 10 },
      tickformat: ".0%",
      rangemode: "tozero",
      gridcolor: "#eee",
    },
    xaxis: {
      tickangle: panels.length > 2 ? -15 : 0,
      automargin: true,
    },
    legend: { orientation: "h", y: -0.25, x: 0.5, xanchor: "center" },
    plot_bgcolor: "#fff",
    paper_bgcolor: "#fff",
    font: { family: "-apple-system, BlinkMacSystemFont, 'Segoe UI', sans-serif", size: 13 },
    height: 380,
    boxmode: "group",
  };

  Plotly.react(containerId, traces, layout, { displayModeBar: false, responsive: true });

  // numeric readout below the plot
  const rows = panels.map((panel) => {
    const a = panel.data.alert;
    const u = panel.data.usual;
    const rr = panel.data.rr == null
      ? "—"
      : `${panel.data.rr.toFixed(2)} (${panel.data.rr_lo.toFixed(2)}–${panel.data.rr_hi.toFixed(2)})`;
    return `
      <tr>
        <td>${escapeHtml(panel.label)}</td>
        <td class="col-num">${a.events}/${a.n} (${fmtPct(a.prop)})</td>
        <td class="col-num">${u.events}/${u.n} (${fmtPct(u.prop)})</td>
        <td class="col-num">${rr}</td>
        <td class="col-num">${fmtP(panel.data.cmh_p)}</td>
      </tr>`;
  }).join("");
  document.getElementById(readoutId).innerHTML = `
    <table class="data-table compact">
      <thead>
        <tr>
          <th>Subgroup</th>
          <th class="col-num">Alert</th>
          <th class="col-num">Usual care</th>
          <th class="col-num">Relative risk (95% CI)</th>
          <th class="col-num">CMH p</th>
        </tr>
      </thead>
      <tbody>${rows}</tbody>
    </table>`;
}

function initFigure(figKey, selectId, plotId, readoutId, captionId, footnoteId, nSummaryId) {
  const fig = state[figKey];
  document.getElementById(captionId).textContent = fig.title;
  document.getElementById(footnoteId).textContent = fig.footnote;

  const select = document.getElementById(selectId);
  populateGroupedSelect(select, fig.subgroups);
  select.addEventListener("change", () => renderFigure(plotId, readoutId, nSummaryId, fig, select.value));
  renderFigure(plotId, readoutId, nSummaryId, fig, "overall");
}

// ---------- Fig 3 (subgroup forest plot) ----------

function fmtPInline(p) {
  if (p == null) return "—";
  if (p < 0.001) return p.toExponential(1).replace("e-", "×10⁻").replace(/×10⁻(\d+)/, (_, d) => "×10⁻" + d);
  if (p < 0.01) return p.toFixed(3);
  return p.toFixed(2);
}

function renderFig3() {
  const fig = state.fig3;
  document.getElementById("fig3-caption").textContent = fig.title;
  document.getElementById("fig3-footnote").textContent = fig.footnote;

  // Build row list: each subgroup contributes a header row + one row per level.
  // We use numeric y-positions (top to bottom, decreasing) so we can place
  // headers and levels freely. The y-axis will use tickvals/ticktext.
  const rows = [];
  let y = 0;
  for (const sg of fig.subgroups) {
    rows.push({ kind: "header", y: y--, name: sg.name, p: sg.interaction_p });
    for (const r of sg.rows) {
      rows.push({ kind: "level", y: y--, ...r });
    }
    y--; // blank spacer
  }

  const levels = rows.filter((r) => r.kind === "level");
  const xs = levels.map((r) => r.rr);
  const ys = levels.map((r) => r.y);
  const errMinus = levels.map((r) => r.rr - r.rr_lo);
  const errPlus = levels.map((r) => r.rr_hi - r.rr);
  const hoverText = levels.map((r) => {
    const a = `${r.events_alert}/${r.n_alert} (${fmtPct(r.p_alert)})`;
    const u = `${r.events_usual}/${r.n_usual} (${fmtPct(r.p_usual)})`;
    return `<b>${r.label}</b><br>Alert: ${a}<br>Usual: ${u}<br>RR ${r.rr.toFixed(2)} (${r.rr_lo.toFixed(2)}–${r.rr_hi.toFixed(2)})`;
  });

  const trace = {
    type: "scatter",
    mode: "markers",
    x: xs,
    y: ys,
    marker: { size: 10, color: ARM_COLORS.alert, symbol: "square" },
    error_x: {
      type: "data", symmetric: false,
      array: errPlus, arrayminus: errMinus,
      color: ARM_COLORS.alert, thickness: 1.5, width: 6,
    },
    text: hoverText,
    hovertemplate: "%{text}<extra></extra>",
    showlegend: false,
  };

  // y-axis tick labels:
  //   header rows   → "<b>Subgroup</b>   interaction p = X"
  //   level rows    → "   Level   ev_a/n_a vs ev_u/n_u"
  const tickvals = rows.map((r) => r.y);
  const ticktext = rows.map((r) => {
    if (r.kind === "header") {
      return `<b>${r.name}</b>   interaction p = ${fmtPInline(r.p)}`;
    }
    const counts = `${r.events_alert}/${r.n_alert} vs ${r.events_usual}/${r.n_usual}`;
    return `   ${r.label}   <span style="color:#6b7280">${counts}</span>`;
  });

  // Right-side annotation per level: RR (lo–hi). Header for the column above.
  const annotations = levels.map((r) => ({
    xref: "paper", x: 1.02, xanchor: "left",
    yref: "y", y: r.y, yanchor: "middle",
    text: `${r.rr.toFixed(2)} (${r.rr_lo.toFixed(2)}–${r.rr_hi.toFixed(2)})`,
    showarrow: false,
    font: { family: "-apple-system, sans-serif", size: 11, color: "#374151" },
  }));
  annotations.push({
    xref: "paper", x: 1.02, xanchor: "left",
    yref: "y", y: 0.6, yanchor: "bottom",
    text: "<b>RR (95% CI)</b>",
    showarrow: false,
    font: { family: "-apple-system, sans-serif", size: 10.5, color: "#374151" },
  });

  const layout = {
    margin: { l: 340, r: 140, t: 22, b: 55 },
    xaxis: {
      type: "log",
      title: { text: "Relative risk (alert vs. usual care, log scale)", standoff: 10 },
      range: [Math.log10(0.5), Math.log10(2.0)],
      tickvals: [0.5, 0.7, 1, 1.5, 2],
      ticktext: ["0.5", "0.7", "1.0", "1.5", "2.0"],
      gridcolor: "#eee",
      zeroline: false,
    },
    yaxis: {
      tickvals,
      ticktext,
      showgrid: false,
      zeroline: false,
      range: [Math.min(...tickvals) - 1, Math.max(...tickvals) + 1],
      tickfont: { family: "-apple-system, sans-serif", size: 12 },
      automargin: false,
    },
    shapes: [
      {
        type: "line",
        xref: "x", x0: 1, x1: 1,
        yref: "paper", y0: 0, y1: 1,
        line: { color: "#9ca3af", width: 1, dash: "dash" },
      },
    ],
    annotations,
    plot_bgcolor: "#fff",
    paper_bgcolor: "#fff",
    font: { family: "-apple-system, BlinkMacSystemFont, 'Segoe UI', sans-serif", size: 13 },
    height: Math.max(320, 28 * rows.length + 80),
  };

  Plotly.react("fig3-plot", [trace], layout, { displayModeBar: false, responsive: true });
}

// ---------- Table 3 (safety outcomes) ----------

function renderTable3() {
  const t3 = state.table3;
  document.getElementById("table3-footnote").textContent = t3.footnote;

  const parts = [];
  parts.push(`
    <table class="data-table">
      <thead>
        <tr>
          <th>Outcome</th>
          <th class="col-num">Alert</th>
          <th class="col-num">Usual care</th>
          <th class="col-num">Difference (95% CI)</th>
        </tr>
      </thead>
      <tbody>`);
  for (const panel of t3.panels) {
    parts.push(`
      <tr class="section-row">
        <td colspan="4">${escapeHtml(panel.label)} · Alert n = ${panel.n_alert.toLocaleString()}, Usual care n = ${panel.n_usual.toLocaleString()}</td>
      </tr>`);
    for (const r of panel.rows) {
      parts.push(`
        <tr>
          <td>${escapeHtml(r.label)}</td>
          <td class="col-num">${escapeHtml(r.alert)}</td>
          <td class="col-num">${escapeHtml(r.usual)}</td>
          <td class="col-num">${escapeHtml(r.diff)}</td>
        </tr>`);
    }
  }
  parts.push(`</tbody></table>`);
  document.getElementById("table3-container").innerHTML = parts.join("");
}

// ---------- Table 2 ----------

function renderTable2(key) {
  const sub = state.table2.subgroups[key];
  document.getElementById("table2-n").textContent =
    `Alert: n = ${sub.n_alert.toLocaleString()} · Usual care: n = ${sub.n_usual.toLocaleString()}`;

  const rows = sub.rows.map((row) => `
    <tr>
      <td>${escapeHtml(row.label)}</td>
      <td class="col-num">${escapeHtml(row.alert)}</td>
      <td class="col-num">${escapeHtml(row.usual)}</td>
      <td class="col-num">${escapeHtml(row.rr)}</td>
      <td class="col-num">${escapeHtml(row.p)}</td>
    </tr>`).join("");

  document.getElementById("table2-container").innerHTML = `
    <table class="data-table">
      <thead>
        <tr>
          <th>Outcome</th>
          <th class="col-num">Alert (N = ${sub.n_alert.toLocaleString()})</th>
          <th class="col-num">Usual care (N = ${sub.n_usual.toLocaleString()})</th>
          <th class="col-num">Relative risk (95% CI)</th>
          <th class="col-num">CMH p</th>
        </tr>
      </thead>
      <tbody>${rows}</tbody>
    </table>`;
}

function initTable2() {
  const select = document.getElementById("table2-subgroup");
  populateGroupedSelect(select, state.table2.subgroups);
  select.addEventListener("change", () => renderTable2(select.value));
  document.getElementById("table2-footnote").textContent = state.table2.footnote;
  renderTable2("overall");
}

// ---------- boot ----------

async function boot() {
  const [t1, f1, f2, f3, t2, t3] = await Promise.all([
    fetchJSON("data/table1.json"),
    fetchJSON("data/fig1.json"),
    fetchJSON("data/fig2.json"),
    fetchJSON("data/fig3.json"),
    fetchJSON("data/table2.json"),
    fetchJSON("data/table3.json"),
  ]);
  state.table1 = t1;
  state.fig1 = f1;
  state.fig2 = f2;
  state.fig3 = f3;
  state.table2 = t2;
  state.table3 = t3;

  initTable1();
  initFigure("fig1", "fig1-view", "fig1-plot", "fig1-readout", "fig1-caption", "fig1-footnote", "fig1-n");
  initFigure("fig2", "fig2-view", "fig2-plot", "fig2-readout", "fig2-caption", "fig2-footnote", "fig2-n");
  initTable2();
  renderFig3();
  renderTable3();
}

boot().catch((err) => {
  console.error(err);
  for (const id of ["table1-container", "fig1-plot", "fig2-plot", "fig3-plot", "table2-container", "table3-container"]) {
    const el = document.getElementById(id);
    if (el) el.textContent = "Failed to load data. If viewing locally, open via a local web server (fetch() requires http://, not file://).";
  }
});
