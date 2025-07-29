let currentAnalysis = null;
let analysisResults = null;

function getElement(id) {
    const element = document.getElementById(id);
    if (!element) {
        console.warn(`Element with id '${id}' not found`);
    }
    return element;
}

async function submitDesign() {
    const sequenceInput = getElement('sequenceInput');
    if (!sequenceInput) {
        showNotification('Sequence input element not found', 'error');
        return;
    }

    const sequence = sequenceInput.value.trim();
    if (!sequence) {
        showNotification('Please enter a sequence', 'error');
        return;
    }

    // Get parameters with null checks
    const parameters = {
        minLength: parseInt(getElement('minLength')?.value || '18'),
        maxLength: parseInt(getElement('maxLength')?.value || '25'),
        minTm: parseFloat(getElement('minTm')?.value || '57.0'),
        maxTm: parseFloat(getElement('maxTm')?.value || '63.0'),
        minGC: parseFloat(getElement('minGC')?.value || '40.0'),
        maxGC: parseFloat(getElement('maxGC')?.value || '60.0'),
        minProduct: parseInt(getElement('minProduct')?.value || '100'),
        maxProduct: parseInt(getElement('maxProduct')?.value || '500')
    };

    // Show progress section
    const inputSection = getElement('input-section');
    const progressSection = getElement('progress-section');
    if (inputSection) inputSection.style.display = 'none';
    if (progressSection) progressSection.style.display = 'block';
    updateStepIndicator(2);

    try {
        // Submit design request
        const response = await fetch('/api/design', {
            method: 'POST',
            headers: {
                'Content-Type': 'application/json'
            },
            body: JSON.stringify({
                sequence: sequence,
                parameters: parameters
            })
        });

        const data = await response.json();
        
        if (!response.ok) {
            throw new Error(data.error || 'Failed to submit design request');
        }

        // Start polling for results
        pollResults(data.job_id);

    } catch (error) {
        showNotification(error.message, 'error');
        resetForm();
    }
}

async function pollResults(jobId) {
    try {
        const response = await fetch(`/api/status/${jobId}`);
        const data = await response.json();

        if (!response.ok) {
            throw new Error(data.error || 'Failed to get job status');
        }

        // Update progress based on status
        if (data.status === 'processing') {
            updateProgressStep(data.progress || 1);
        }

        if (data.status === 'completed') {
            // Process and display results
            displayResults(data.results);
        } else if (data.status === 'failed') {
            throw new Error(data.error_message || 'Primer design failed');
        } else {
            // Continue polling with a longer interval
            setTimeout(() => pollResults(jobId), 5000);
        }
    } catch (error) {
        showNotification(error.message, 'error');
        resetForm();
    }
}

function updateProgressStep(step) {
    const steps = document.querySelectorAll('.progress-step');
    if (steps.length === 0) return;
    
    steps.forEach((el, index) => {
        if (index < step) {
            el.style.display = 'block';
        } else {
            el.style.display = 'none';
        }
    });
}

function displayResults(results) {
    // Store results in window object for access by action buttons
    window.currentResults = results;
    
    if (!results || !results.primers || results.primers.length === 0) {
        showNotification('No primers could be designed', 'error');
        resetForm();
        return;
    }

    const primers = results.primers;
    const resultsContainer = getElement('resultsContainer');
    if (!resultsContainer) return;

    // Update summary stats
    const totalPrimers = getElement('totalPrimers');
    const highQuality = getElement('highQuality');
    const mediumQuality = getElement('mediumQuality');
    const averageTm = getElement('averageTm');

    // Count quality levels based on overall_score
    const qualityCounts = {
        high: primers.filter(p => p.quality_scores?.overall_score >= 80).length,
        medium: primers.filter(p => p.quality_scores?.overall_score >= 60 && p.quality_scores?.overall_score < 80).length,
        low: primers.filter(p => p.quality_scores?.overall_score < 60).length
    };

    if (totalPrimers) totalPrimers.textContent = primers.length;
    if (highQuality) highQuality.textContent = qualityCounts.high;
    if (mediumQuality) mediumQuality.textContent = qualityCounts.medium;
    if (averageTm) {
        const avgTm = primers.reduce((sum, p) => sum + p.forward_tm, 0) / primers.length;
        averageTm.textContent = `${avgTm.toFixed(1)}°C`;
    }

    // Create results table
    const table = document.createElement('table');
    table.className = 'table table-striped table-hover';
    table.innerHTML = `
        <thead>
            <tr>
                <th>Pair</th>
                <th>Forward Primer</th>
                <th>Reverse Primer</th>
                <th>Product Size</th>
                <th>TM</th>
                <th>GC%</th>
                <th>Quality</th>
                <th>Actions</th>
            </tr>
        </thead>
        <tbody></tbody>
    `;

    const tbody = table.querySelector('tbody');
    primers.forEach((primer, index) => {
        // Get quality scores
        const qualityScores = primer.quality_scores || {};
        const overallScore = qualityScores.overall_score || 0;
        
        // Determine quality level and color
        let qualityLevel = 'low';
        let qualityColor = 'danger';
        if (overallScore >= 80) {
            qualityLevel = 'high';
            qualityColor = 'success';
        } else if (overallScore >= 60) {
            qualityLevel = 'medium';
            qualityColor = 'warning';
        }

        // Get individual scores
        const lengthScore = qualityScores.length_score || 0;
        const gcScore = qualityScores.gc_score || 0;
        const tmScore = qualityScores.tm_score || 0;
        const complexityScore = qualityScores.complexity_score || 0;
        const terminalScore = qualityScores.terminal_score || 0;

        const row = document.createElement('tr');
        row.className = 'primer-row';
        row.setAttribute('data-forward', primer.forward_primer);
        row.setAttribute('data-reverse', primer.reverse_primer);
        row.innerHTML = `
            <td>${index + 1}</td>
            <td>
                <div class="primer-sequence">${primer.forward_primer}</div>
                <small class="text-muted">${primer.forward_position}</small>
            </td>
            <td>
                <div class="primer-sequence">${primer.reverse_primer}</div>
                <small class="text-muted">${primer.reverse_position}</small>
            </td>
            <td>${primer.product_size} bp</td>
            <td>${primer.forward_tm.toFixed(1)}°C</td>
            <td>${primer.gc_content}%</td>
            <td>
                <span class="badge bg-${qualityColor}">
                    ${qualityLevel.toUpperCase()}
                </span>
                <small class="d-block text-muted">Overall: ${overallScore.toFixed(1)}</small>
                <small class="d-block text-muted">Length: ${lengthScore.toFixed(1)}</small>
                <small class="d-block text-muted">GC: ${gcScore.toFixed(1)}</small>
                <small class="d-block text-muted">TM: ${tmScore.toFixed(1)}</small>
                <small class="d-block text-muted">Complexity: ${complexityScore.toFixed(1)}</small>
                <small class="d-block text-muted">Terminal: ${terminalScore.toFixed(1)}</small>
            </td>
            <td>
                <button class="btn btn-sm btn-outline-primary" onclick="showDetails(${index})">
                    <i class="fas fa-info-circle"></i>
                </button>
                <button class="btn btn-sm btn-outline-success" onclick="copyPrimers(${index})">
                    <i class="fas fa-copy"></i>
                </button>
            </td>
        `;
        tbody.appendChild(row);
    });

    // Clear existing content
    resultsContainer.innerHTML = '';
    resultsContainer.appendChild(table);

    // Clear and create visualizations
    const bindingSitesViz = getElement('bindingSitesViz');
    const meltingProfileViz = getElement('meltingProfileViz');
    const dimerViz = getElement('dimerViz');
    const qualityChart = getElement('qualityChart');

    if (bindingSitesViz) {
        bindingSitesViz.innerHTML = '';
        createBindingSitesViz(primers);
    }
    if (meltingProfileViz) {
        meltingProfileViz.innerHTML = '';
        createMeltingProfileViz(primers);
    }
    if (dimerViz) {
        dimerViz.innerHTML = '';
        createDimerVisualizations(primers);
    }
    if (qualityChart) {
        qualityChart.innerHTML = '';
        createQualityChart(primers);
    }

    // Show results section with animation
    const progressSection = getElement('progress-section');
    const resultsSection = getElement('results-section');
    if (progressSection) {
        progressSection.classList.remove('animate__fadeIn');
        progressSection.classList.add('animate__fadeOut');
        setTimeout(() => {
            progressSection.style.display = 'none';
            progressSection.classList.remove('animate__fadeOut');
        }, 500);
    }
    if (resultsSection) {
        resultsSection.style.display = 'block';
        resultsSection.classList.remove('animate__fadeOutUp');
        resultsSection.classList.add('animate__fadeInUp');
    }
    updateStepIndicator(4);

    // Update cost estimate
    updateCostEstimate();
}

function createBindingSitesViz(primers) {
    const container = getElement('bindingSitesViz');
    if (!container) return;

    // Create a simple visualization of primer binding sites
    const sequenceLength = Math.max(...primers.map(p => Math.max(p.forward_position, p.reverse_position)));
    const viz = document.createElement('div');
    viz.style.position = 'relative';
    viz.style.height = '200px';
    viz.style.border = '1px solid #ddd';
    viz.style.borderRadius = '4px';
    viz.style.padding = '10px';
    viz.style.overflow = 'auto';

    // Create sequence line
    const sequenceLine = document.createElement('div');
    sequenceLine.style.position = 'relative';
    sequenceLine.style.height = '2px';
    sequenceLine.style.background = '#666';
    sequenceLine.style.margin = '50px 0';
    viz.appendChild(sequenceLine);

    // Add primer binding sites
    primers.forEach((primer, index) => {
        // Forward primer
        const forward = document.createElement('div');
        forward.style.position = 'absolute';
        forward.style.left = `${(primer.forward_position / sequenceLength) * 100}%`;
        forward.style.top = '20px';
        forward.style.width = '4px';
        forward.style.height = '20px';
        forward.style.background = '#007bff';
        forward.title = `Forward Primer ${index + 1}: ${primer.forward_primer}`;
        viz.appendChild(forward);

        // Reverse primer
        const reverse = document.createElement('div');
        reverse.style.position = 'absolute';
        reverse.style.left = `${(primer.reverse_position / sequenceLength) * 100}%`;
        reverse.style.top = '60px';
        reverse.style.width = '4px';
        reverse.style.height = '20px';
        reverse.style.background = '#dc3545';
        reverse.title = `Reverse Primer ${index + 1}: ${primer.reverse_primer}`;
        viz.appendChild(reverse);

        // Connect primers with a line
        const line = document.createElement('div');
        line.style.position = 'absolute';
        line.style.left = `${(primer.forward_position / sequenceLength) * 100}%`;
        line.style.top = '30px';
        line.style.width = `${((primer.reverse_position - primer.forward_position) / sequenceLength) * 100}%`;
        line.style.height = '2px';
        line.style.background = '#28a745';
        line.title = `Product Size: ${primer.product_size} bp`;
        viz.appendChild(line);
    });

    container.innerHTML = '';
    container.appendChild(viz);
}

function createMeltingProfileViz(primers) {
    const container = getElement('meltingProfileViz');
    if (!container) return;

    // Create a simple bar chart of melting temperatures
    const viz = document.createElement('div');
    viz.style.position = 'relative';
    viz.style.height = '200px';
    viz.style.border = '1px solid #ddd';
    viz.style.borderRadius = '4px';
    viz.style.padding = '10px';

    // Calculate min and max temperatures
    const minTemp = Math.min(...primers.map(p => p.forward_tm));
    const maxTemp = Math.max(...primers.map(p => p.forward_tm));
    const range = maxTemp - minTemp;

    // Create bars for each primer pair
    primers.forEach((primer, index) => {
        const bar = document.createElement('div');
        bar.style.position = 'absolute';
        bar.style.left = `${(index / primers.length) * 100}%`;
        bar.style.bottom = '0';
        bar.style.width = `${(1 / primers.length) * 90}%`;
        bar.style.height = `${((primer.forward_tm - minTemp) / range) * 80}%`;
        bar.style.background = '#007bff';
        bar.style.margin = '0 2px';
        bar.title = `Pair ${index + 1}: ${primer.forward_tm.toFixed(1)}°C`;
        viz.appendChild(bar);
    });

    // Add temperature scale
    const scale = document.createElement('div');
    scale.style.position = 'absolute';
    scale.style.left = '0';
    scale.style.right = '0';
    scale.style.bottom = '-20px';
    scale.style.textAlign = 'center';
    scale.style.fontSize = '12px';
    scale.style.color = '#666';
    scale.textContent = `${minTemp.toFixed(1)}°C - ${maxTemp.toFixed(1)}°C`;
    viz.appendChild(scale);

    container.innerHTML = '';
    container.appendChild(viz);
}

function createDimerVisualizations(primers) {
    const container = getElement('dimerViz');
    if (!container) return;

    // Create a grid of primer pair interactions
    const viz = document.createElement('div');
    viz.style.display = 'grid';
    viz.style.gridTemplateColumns = `repeat(${primers.length}, 1fr)`;
    viz.style.gap = '5px';
    viz.style.padding = '10px';

    // Create header row
    const headerRow = document.createElement('div');
    headerRow.style.gridColumn = '1 / -1';
    headerRow.style.display = 'grid';
    headerRow.style.gridTemplateColumns = `repeat(${primers.length}, 1fr)`;
    headerRow.style.gap = '5px';
    headerRow.style.marginBottom = '5px';

    // Add column headers
    for (let i = 0; i < primers.length; i++) {
        const header = document.createElement('div');
        header.style.textAlign = 'center';
        header.style.fontSize = '12px';
        header.style.fontWeight = 'bold';
        header.textContent = `Pair ${i + 1}`;
        headerRow.appendChild(header);
    }
    viz.appendChild(headerRow);

    // Create interaction matrix
    for (let i = 0; i < primers.length; i++) {
        for (let j = 0; j < primers.length; j++) {
            const cell = document.createElement('div');
            cell.style.padding = '5px';
            cell.style.textAlign = 'center';
            cell.style.border = '1px solid #ddd';
            cell.style.borderRadius = '4px';
            
            if (i === j) {
                cell.style.background = '#f8f9fa';
                cell.textContent = '-';
            } else {
                // Calculate potential dimer score (simplified)
                const score = Math.random() * 100; // Replace with actual dimer calculation
                cell.style.background = score > 70 ? '#dc3545' : score > 40 ? '#ffc107' : '#28a745';
                cell.textContent = `${score.toFixed(1)}%`;
                cell.title = `Dimer score between Pair ${i + 1} and Pair ${j + 1}`;
            }
            viz.appendChild(cell);
        }
    }

    container.innerHTML = '';
    container.appendChild(viz);
}

function createQualityChart(primers) {
    const container = getElement('qualityChart');
    if (!container) return;

    // Create a simple pie chart of quality distribution
    const viz = document.createElement('div');
    viz.style.position = 'relative';
    viz.style.height = '300px';
    viz.style.marginBottom = '60px'; // Add margin to prevent overlap with button
    viz.style.border = '1px solid #ddd';
    viz.style.borderRadius = '4px';
    viz.style.padding = '10px';
    viz.style.backgroundColor = '#fff';

    // Count quality levels based on overall_score
    const qualityCounts = {
        high: primers.filter(p => p.quality_scores?.overall_score >= 80).length,
        medium: primers.filter(p => p.quality_scores?.overall_score >= 60 && p.quality_scores?.overall_score < 80).length,
        low: primers.filter(p => p.quality_scores?.overall_score < 60).length
    };

    // Create pie chart
    const total = primers.length;
    let startAngle = 0;
    const colors = {
        high: '#28a745',
        medium: '#ffc107',
        low: '#dc3545'
    };

    // Create SVG element
    const svg = document.createElementNS('http://www.w3.org/2000/svg', 'svg');
    svg.setAttribute('width', '100%');
    svg.setAttribute('height', '100%');
    svg.style.position = 'absolute';
    svg.style.top = '0';
    svg.style.left = '0';

    // Calculate center and radius
    const width = viz.clientWidth;
    const height = viz.clientHeight;
    const centerX = width / 2;
    const centerY = height / 2;
    const radius = Math.min(width, height) * 0.35; // Slightly reduce radius to prevent overlap

    // Draw pie segments
    Object.entries(qualityCounts).forEach(([level, count]) => {
        if (count === 0) return;

        const percentage = (count / total) * 100;
        const angle = (percentage / 100) * 360;
        const endAngle = startAngle + angle;

        // Convert angles to radians
        const startRad = (startAngle - 90) * Math.PI / 180;
        const endRad = (endAngle - 90) * Math.PI / 180;

        // Calculate points for the arc
        const x1 = centerX + radius * Math.cos(startRad);
        const y1 = centerY + radius * Math.sin(startRad);
        const x2 = centerX + radius * Math.cos(endRad);
        const y2 = centerY + radius * Math.sin(endRad);

        // Create path for the pie segment
        const path = document.createElementNS('http://www.w3.org/2000/svg', 'path');
        const largeArcFlag = angle > 180 ? 1 : 0;
        const pathData = [
            `M ${centerX},${centerY}`,
            `L ${x1},${y1}`,
            `A ${radius},${radius} 0 ${largeArcFlag} 1 ${x2},${y2}`,
            'Z'
        ].join(' ');

        path.setAttribute('d', pathData);
        path.setAttribute('fill', colors[level]);
        path.setAttribute('stroke', '#fff');
        path.setAttribute('stroke-width', '2');

        // Add hover effect
        path.style.transition = 'opacity 0.3s';
        path.addEventListener('mouseover', () => {
            path.style.opacity = '0.8';
        });
        path.addEventListener('mouseout', () => {
            path.style.opacity = '1';
        });

        svg.appendChild(path);

        // Add label
        const labelAngle = (startAngle + endAngle) / 2;
        const labelRad = (labelAngle - 90) * Math.PI / 180;
        const labelX = centerX + (radius * 0.6) * Math.cos(labelRad);
        const labelY = centerY + (radius * 0.6) * Math.sin(labelRad);

        const text = document.createElementNS('http://www.w3.org/2000/svg', 'text');
        text.setAttribute('x', labelX);
        text.setAttribute('y', labelY);
        text.setAttribute('text-anchor', 'middle');
        text.setAttribute('dominant-baseline', 'middle');
        text.setAttribute('fill', '#fff');
        text.setAttribute('font-weight', 'bold');
        text.textContent = `${Math.round(percentage)}%`;
        svg.appendChild(text);

        startAngle = endAngle;
    });

    // Add legend
    const legend = document.createElement('div');
    legend.style.position = 'absolute';
    legend.style.bottom = '10px';
    legend.style.left = '50%';
    legend.style.transform = 'translateX(-50%)';
    legend.style.display = 'flex';
    legend.style.gap = '20px';
    legend.style.backgroundColor = '#fff';
    legend.style.padding = '5px 10px';
    legend.style.borderRadius = '4px';
    legend.style.boxShadow = '0 1px 3px rgba(0,0,0,0.1)';

    Object.entries(qualityCounts).forEach(([level, count]) => {
        const legendItem = document.createElement('div');
        legendItem.style.display = 'flex';
        legendItem.style.alignItems = 'center';
        legendItem.style.gap = '5px';

        const colorBox = document.createElement('div');
        colorBox.style.width = '12px';
        colorBox.style.height = '12px';
        colorBox.style.backgroundColor = colors[level];
        colorBox.style.borderRadius = '2px';

        const label = document.createElement('span');
        label.textContent = `${level.charAt(0).toUpperCase() + level.slice(1)} (${count})`;

        legendItem.appendChild(colorBox);
        legendItem.appendChild(label);
        legend.appendChild(legendItem);
    });

    viz.appendChild(svg);
    viz.appendChild(legend);
    container.appendChild(viz);
}

function downloadViz(type) {
    // Implement visualization download functionality
    showNotification('Download functionality will be implemented', 'info');
}

function showDetails(index) {
    const primers = window.currentResults?.primers;
    if (!primers || !primers[index]) return;

    const primer = primers[index];
    const qualityScores = primer.quality_scores || {};

    // Create modal content
    const modalContent = `
        <div class="modal fade" id="primerDetailsModal" tabindex="-1">
            <div class="modal-dialog modal-lg">
                <div class="modal-content">
                    <div class="modal-header">
                        <h5 class="modal-title">Primer Pair Details</h5>
                        <button type="button" class="btn-close" data-bs-dismiss="modal"></button>
                    </div>
                    <div class="modal-body">
                        <div class="row mb-3">
                            <div class="col-md-6">
                                <h6>Forward Primer</h6>
                                <p><strong>Sequence:</strong> ${primer.forward_primer}</p>
                                <p><strong>Position:</strong> ${primer.forward_position}</p>
                                <p><strong>Tm:</strong> ${primer.forward_tm.toFixed(1)}°C</p>
                    </div>
                    <div class="col-md-6">
                                <h6>Reverse Primer</h6>
                                <p><strong>Sequence:</strong> ${primer.reverse_primer}</p>
                                <p><strong>Position:</strong> ${primer.reverse_position}</p>
                                <p><strong>Tm:</strong> ${primer.reverse_tm.toFixed(1)}°C</p>
                            </div>
                        </div>
                        <div class="row mb-3">
                            <div class="col-12">
                                <h6>Product Information</h6>
                                <p><strong>Size:</strong> ${primer.product_size} bp</p>
                        <p><strong>GC Content:</strong> ${primer.gc_content}%</p>
                            </div>
                        </div>
                        <div class="row">
                            <div class="col-12">
                                <h6>Quality Scores</h6>
                                <div class="table-responsive">
                                    <table class="table table-sm">
                                        <thead>
                                            <tr>
                                                <th>Metric</th>
                                                <th>Score</th>
                                                <th>Status</th>
                                            </tr>
                                        </thead>
                                        <tbody>
                                            <tr>
                                                <td>Overall Quality</td>
                                                <td>${qualityScores.overall_score?.toFixed(1) || 'N/A'}</td>
                                                <td>
                                                    <span class="badge bg-${qualityScores.overall_score >= 80 ? 'success' : 
                                                                         qualityScores.overall_score >= 60 ? 'warning' : 
                                                                         'danger'}">
                                                        ${qualityScores.overall_score >= 80 ? 'High' : 
                                                          qualityScores.overall_score >= 60 ? 'Medium' : 
                                                          'Low'}
                                                    </span>
                                                </td>
                                            </tr>
                                            <tr>
                                                <td>Length Score</td>
                                                <td>${qualityScores.length_score?.toFixed(1) || 'N/A'}</td>
                                                <td>
                                                    <span class="badge bg-${qualityScores.length_score >= 80 ? 'success' : 
                                                                         qualityScores.length_score >= 60 ? 'warning' : 
                                                                         'danger'}">
                                                        ${qualityScores.length_score >= 80 ? 'Good' : 
                                                          qualityScores.length_score >= 60 ? 'Fair' : 
                                                          'Poor'}
                                                    </span>
                                                </td>
                                            </tr>
                                            <tr>
                                                <td>GC Content Score</td>
                                                <td>${qualityScores.gc_score?.toFixed(1) || 'N/A'}</td>
                                                <td>
                                                    <span class="badge bg-${qualityScores.gc_score >= 80 ? 'success' : 
                                                                         qualityScores.gc_score >= 60 ? 'warning' : 
                                                                         'danger'}">
                                                        ${qualityScores.gc_score >= 80 ? 'Good' : 
                                                          qualityScores.gc_score >= 60 ? 'Fair' : 
                                                          'Poor'}
                                                    </span>
                                                </td>
                                            </tr>
                                            <tr>
                                                <td>Melting Temperature Score</td>
                                                <td>${qualityScores.tm_score?.toFixed(1) || 'N/A'}</td>
                                                <td>
                                                    <span class="badge bg-${qualityScores.tm_score >= 80 ? 'success' : 
                                                                         qualityScores.tm_score >= 60 ? 'warning' : 
                                                                         'danger'}">
                                                        ${qualityScores.tm_score >= 80 ? 'Good' : 
                                                          qualityScores.tm_score >= 60 ? 'Fair' : 
                                                          'Poor'}
                                                    </span>
                                                </td>
                                            </tr>
                                            <tr>
                                                <td>Complexity Score</td>
                                                <td>${qualityScores.complexity_score?.toFixed(1) || 'N/A'}</td>
                                                <td>
                                                    <span class="badge bg-${qualityScores.complexity_score >= 80 ? 'success' : 
                                                                         qualityScores.complexity_score >= 60 ? 'warning' : 
                                                                         'danger'}">
                                                        ${qualityScores.complexity_score >= 80 ? 'Good' : 
                                                          qualityScores.complexity_score >= 60 ? 'Fair' : 
                                                          'Poor'}
                                                    </span>
                                                </td>
                                            </tr>
                                            <tr>
                                                <td>Terminal Score</td>
                                                <td>${qualityScores.terminal_score?.toFixed(1) || 'N/A'}</td>
                                                <td>
                                                    <span class="badge bg-${qualityScores.terminal_score >= 80 ? 'success' : 
                                                                         qualityScores.terminal_score >= 60 ? 'warning' : 
                                                                         'danger'}">
                                                        ${qualityScores.terminal_score >= 80 ? 'Good' : 
                                                          qualityScores.terminal_score >= 60 ? 'Fair' : 
                                                          'Poor'}
                                                    </span>
                                                </td>
                                            </tr>
                                        </tbody>
                                    </table>
                                </div>
                            </div>
                        </div>
                    </div>
                    <div class="modal-footer">
                        <button type="button" class="btn btn-secondary" data-bs-dismiss="modal">Close</button>
                    </div>
                </div>
            </div>
            </div>
        `;

    // Remove existing modal if any
    const existingModal = document.getElementById('primerDetailsModal');
    if (existingModal) {
        existingModal.remove();
    }

    // Add new modal to document
    document.body.insertAdjacentHTML('beforeend', modalContent);

    // Initialize and show modal
    const modal = new bootstrap.Modal(document.getElementById('primerDetailsModal'));
    modal.show();
}

function copyPrimers(index) {
    const primers = window.currentResults?.primers;
    if (!primers || !primers[index]) return;

    const primer = primers[index];
    const textToCopy = `Forward Primer: ${primer.forward_primer}
Reverse Primer: ${primer.reverse_primer}
Product Size: ${primer.product_size} bp
Tm: ${primer.forward_tm.toFixed(1)}°C (F), ${primer.reverse_tm.toFixed(1)}°C (R)
GC Content: ${primer.gc_content}%
Quality Score: ${primer.quality_scores?.overall_score?.toFixed(1) || 'N/A'}`;

    // Copy to clipboard
    navigator.clipboard.writeText(textToCopy).then(() => {
        showNotification('Primer information copied to clipboard!', 'success');
    }).catch(err => {
        showNotification('Failed to copy primer information', 'error');
        console.error('Failed to copy text: ', err);
    });
}

function showNotification(message, type = 'info', isHtml = false) {
    let notification = document.getElementById('notification');
    if (!notification) {
        notification = document.createElement('div');
        notification.id = 'notification';
        notification.style.position = 'fixed';
        notification.style.top = '30px';
        notification.style.right = '30px';
        notification.style.zIndex = '9999';
        notification.style.minWidth = '280px';
        notification.className = 'alert animate__animated';
        document.body.appendChild(notification);
    }
    notification.className = `alert alert-${type} animate__animated animate__fadeInDown`;
    notification.innerHTML = isHtml ? message : `<span>${message}</span>`;
    notification.style.display = 'block';
    setTimeout(() => {
        notification.classList.remove('animate__fadeInDown');
        notification.classList.add('animate__fadeOutUp');
        setTimeout(() => {
            notification.style.display = 'none';
            notification.classList.remove('animate__fadeOutUp');
        }, 800);
    }, 2500);
}

function updateStepIndicator(step) {
    const steps = document.querySelectorAll('.step');
    if (steps.length === 0) return;

    steps.forEach((el, index) => {
        if (index + 1 < step) {
            el.classList.add('completed');
            el.classList.remove('active');
        } else if (index + 1 === step) {
            el.classList.add('active');
            el.classList.remove('completed');
        } else {
            el.classList.remove('active', 'completed');
        }
    });
}

function resetForm() {
    // Reset form fields
    const sequenceInput = getElement('sequence');
    const minLengthInput = getElement('minLength');
    const maxLengthInput = getElement('maxLength');
    const minTmInput = getElement('minTm');
    const maxTmInput = getElement('maxTm');
    const minGcInput = getElement('minGC');
    const maxGcInput = getElement('maxGC');
    const numReturnInput = getElement('numReturn');

    if (sequenceInput) sequenceInput.value = '';
    if (minLengthInput) minLengthInput.value = '18';
    if (maxLengthInput) maxLengthInput.value = '25';
    if (minTmInput) minTmInput.value = '55';
    if (maxTmInput) maxTmInput.value = '65';
    if (minGcInput) minGcInput.value = '40';
    if (maxGcInput) maxGcInput.value = '60';
    if (numReturnInput) numReturnInput.value = '5';

    // Reset synthesis parameters
    const synthesisScale = getElement('synthesisScale');
    const purificationMethod = getElement('purificationMethod');
    if (synthesisScale) synthesisScale.value = '25nm';
    if (purificationMethod) purificationMethod.value = 'desalt';

    // Clear results
    const resultsContainer = getElement('resultsContainer');
    if (resultsContainer) resultsContainer.innerHTML = '';

    // Reset summary stats
    const totalPrimers = getElement('totalPrimers');
    const highQuality = getElement('highQuality');
    const mediumQuality = getElement('mediumQuality');
    const averageTm = getElement('averageTm');
    if (totalPrimers) totalPrimers.textContent = '0';
    if (highQuality) highQuality.textContent = '0';
    if (mediumQuality) mediumQuality.textContent = '0';
    if (averageTm) averageTm.textContent = '0°C';

    // Clear visualizations
    const bindingSitesViz = getElement('bindingSitesViz');
    const meltingProfileViz = getElement('meltingProfileViz');
    const dimerViz = getElement('dimerViz');
    const qualityChart = getElement('qualityChart');
    if (bindingSitesViz) bindingSitesViz.innerHTML = '';
    if (meltingProfileViz) meltingProfileViz.innerHTML = '';
    if (dimerViz) dimerViz.innerHTML = '';
    if (qualityChart) qualityChart.innerHTML = '';

    // Reset cost estimate
    const costEstimate = getElement('costEstimate');
    if (costEstimate) costEstimate.textContent = '$0.00';

    // Animate hiding results section and showing input section
    const resultsSection = getElement('results-section');
    const inputSection = getElement('input-section');
    if (resultsSection) {
        resultsSection.classList.remove('animate__fadeInUp');
        resultsSection.classList.add('animate__fadeOutDown');
        setTimeout(() => {
            resultsSection.style.display = 'none';
            resultsSection.classList.remove('animate__fadeOutDown');
            if (inputSection) {
                inputSection.style.display = 'block';
                inputSection.classList.add('animate__fadeInLeft');
            }
        }, 500);
    } else if (inputSection) {
        inputSection.style.display = 'block';
        inputSection.classList.add('animate__fadeInLeft');
    }
    
    // Reset step indicator
    updateStepIndicator(1);
    
    // Clear any notifications
    const notificationContainer = getElement('notificationContainer');
    if (notificationContainer) notificationContainer.innerHTML = '';

    // Clear current results
    window.currentResults = null;

    // Focus on sequence input
    if (sequenceInput) sequenceInput.focus();
}

function updateCostEstimate() {
    const synthesisScale = getElement('synthesisScale');
    const purificationMethod = getElement('purificationMethod');
    const costEstimate = getElement('costEstimate');
    
    if (!synthesisScale || !purificationMethod || !costEstimate) return;

    const scale = synthesisScale.value;
    const purification = purificationMethod.value;
    
    // Get all primer rows
    const primerRows = document.querySelectorAll('.primer-row');
    if (!primerRows.length) {
        costEstimate.textContent = '$0.00';
        console.log('[COST DEBUG] No primer rows found.');
        return;
    }

    let totalCost = 0;
    primerRows.forEach(row => {
        const forwardSeq = row.getAttribute('data-forward');
        const reverseSeq = row.getAttribute('data-reverse');
        console.log('[COST DEBUG] forwardSeq:', forwardSeq, 'reverseSeq:', reverseSeq);
        if (forwardSeq) {
            totalCost += calculatePrimerCost(forwardSeq, scale, purification);
        }
        if (reverseSeq) {
            totalCost += calculatePrimerCost(reverseSeq, scale, purification);
        }
    });

    // Format cost with 2 decimal places
    costEstimate.textContent = `$${totalCost.toFixed(2)}`;
    console.log('[COST DEBUG] totalCost:', totalCost);
}

function calculatePrimerCost(sequence, scale, purification) {
    if (!sequence) return 0;

    // Base costs per base pair for different scales
    const baseCosts = {
        '25nm': 0.25,    // $0.25 per base
        '100nm': 0.20,   // $0.20 per base
        '250nm': 0.15,   // $0.15 per base
        '1um': 0.10      // $0.10 per base
    };

    // Additional costs for purification methods
    const purificationCosts = {
        'desalt': 5.00,      // $5.00 per primer
        'standard': 5.00,    // $5.00 per primer (treat 'standard' as 'desalt')
        'hplc': 15.00,       // $15.00 per primer
        'page': 25.00        // $25.00 per primer
    };

    // Calculate base cost
    const baseCost = baseCosts[scale] || baseCosts['25nm'];
    const sequenceLength = sequence.length;
    let cost = sequenceLength * baseCost;

    // Add purification cost
    cost += purificationCosts[purification] || purificationCosts['desalt'];

    // Add minimum order cost
    const minimumOrderCost = 10.00; // $10.00 minimum order cost
    cost = Math.max(cost, minimumOrderCost);

    return cost;
}

function exportResults(format) {
    const results = window.currentResults;
    const devInfo = 'Developed by Manzoor Mohammad Loriya | loryamanjue786@gmail.com | © 2025';
    if (!results || !results.primers || results.primers.length === 0) {
        showNotification('No results to export', 'warning');
        return;
    }
    if (format === 'csv') {
        // Convert results to CSV
        const headers = [
            'Pair', 'Forward Primer', 'Reverse Primer', 'Product Size', 'Forward Tm', 'Reverse Tm', 'GC%',
            'Overall Quality', 'Length Score', 'GC Score', 'Tm Score', 'Complexity Score', 'Terminal Score'
        ];
        const rows = results.primers.map((primer, idx) => [
            idx + 1,
            primer.forward_primer,
            primer.reverse_primer,
            primer.product_size,
            primer.forward_tm,
            primer.reverse_tm,
            primer.gc_content,
            primer.quality_scores?.overall_score ?? '',
            primer.quality_scores?.length_score ?? '',
            primer.quality_scores?.gc_score ?? '',
            primer.quality_scores?.tm_score ?? '',
            primer.quality_scores?.complexity_score ?? '',
            primer.quality_scores?.terminal_score ?? ''
        ]);
        let csv = `# ${devInfo}\n`;
        csv += headers.join(',') + '\n';
        rows.forEach(row => {
            csv += row.map(val => `"${val}"`).join(',') + '\n';
        });
        const blob = new Blob([csv], { type: 'text/csv' });
        const url = URL.createObjectURL(blob);
        const a = document.createElement('a');
        a.href = url;
        a.download = 'primer_results.csv';
        document.body.appendChild(a);
        a.click();
        document.body.removeChild(a);
        URL.revokeObjectURL(url);
        showNotification('Results exported as CSV', 'success');
    } else if (format === 'json') {
        const jsonHeader = `/* ${devInfo} */\n`;
        const json = jsonHeader + JSON.stringify(results, null, 2);
        const blob = new Blob([json], { type: 'application/json' });
        const url = URL.createObjectURL(blob);
        const a = document.createElement('a');
        a.href = url;
        a.download = 'primer_results.json';
        document.body.appendChild(a);
        a.click();
        document.body.removeChild(a);
        URL.revokeObjectURL(url);
        showNotification('Results exported as JSON', 'success');
    } else if (format === 'pdf') {
        if (typeof window.jspdf === 'undefined' && typeof window.jsPDF === 'undefined') {
            showNotification('PDF export requires jsPDF. Please refresh the page.', 'error');
            return;
        }
        // Use jsPDF
        const doc = new (window.jsPDF || window.jspdf.jsPDF)({ orientation: 'landscape' });
        doc.setFontSize(14);
        doc.text('Primer Design Results', 14, 15);
        doc.setFontSize(10);
        doc.text(devInfo, 14, 22);
        const headers = [
            'Pair', 'Forward Primer', 'Reverse Primer', 'Product Size', 'Fwd Tm', 'Rev Tm', 'GC%', 'Quality'
        ];
        const rows = results.primers.map((primer, idx) => [
            idx + 1,
            primer.forward_primer,
            primer.reverse_primer,
            primer.product_size,
            primer.forward_tm,
            primer.reverse_tm,
            primer.gc_content,
            primer.quality_scores?.overall_score ?? ''
        ]);
        // AutoTable if available
        if (doc.autoTable) {
            doc.autoTable({ head: [headers], body: rows, startY: 28, styles: { fontSize: 10 } });
        } else {
            // Fallback: simple text table
            let y = 28;
            doc.setFontSize(10);
            doc.text(headers.join(' | '), 14, y);
            y += 7;
            rows.forEach(row => {
                doc.text(row.join(' | '), 14, y);
                y += 7;
            });
        }
        doc.save('primer_results.pdf');
        showNotification('Results exported as PDF', 'success');
    } else {
        showNotification('Unknown export format', 'error');
    }
}

// Add event listeners when document is loaded
document.addEventListener('DOMContentLoaded', function() {
    const designForm = getElement('designForm');
    if (designForm) {
        designForm.addEventListener('submit', function(e) {
        e.preventDefault();
        submitDesign();
    });
    }

    const synthesisScale = getElement('synthesisScale');
    const purificationMethod = getElement('purificationMethod');
    
    if (synthesisScale) {
        synthesisScale.addEventListener('change', updateCostEstimate);
    }
    if (purificationMethod) {
        purificationMethod.addEventListener('change', updateCostEstimate);
    }

    // Add event listener for the "Design New Primers" button
    const newDesignBtn = getElement('newDesignBtn');
    if (newDesignBtn) {
        newDesignBtn.addEventListener('click', resetForm);
    }

    // Dark mode toggle logic
    const darkModeToggle = document.getElementById('darkModeToggle');
    function setDarkMode(enabled) {
        if (enabled) {
            document.body.classList.add('dark-mode');
            localStorage.setItem('darkMode', '1');
            if (darkModeToggle) darkModeToggle.innerHTML = '<i class="fas fa-sun"></i>';
        } else {
            document.body.classList.remove('dark-mode');
            localStorage.setItem('darkMode', '0');
            if (darkModeToggle) darkModeToggle.innerHTML = '<i class="fas fa-moon"></i>';
        }
    }
    if (darkModeToggle) {
        darkModeToggle.addEventListener('click', function() {
            setDarkMode(!document.body.classList.contains('dark-mode'));
        });
    }
    // On load, set dark mode from localStorage
    if (localStorage.getItem('darkMode') === '1') {
        setDarkMode(true);
    }
}); 