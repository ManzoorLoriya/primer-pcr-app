<<<<<<< HEAD
import plotly.graph_objects as go
from plotly.subplots import make_subplots

def calculate_gc_content(sequence):
    """Calculate GC content percentage of a sequence"""
    gc_count = sequence.count('G') + sequence.count('C')
    return (gc_count / len(sequence)) * 100 if sequence else 0

from .algorithms import PrimerAnalyzer

class VisualizationEngine:
    def __init__(self):
        self.analyzer = PrimerAnalyzer()

    def plot_binding_sites(self, sequence, primers):
        """Visualize primer binding positions on target sequence"""
        fig = go.Figure()
        
        # Add sequence as x-axis
        positions = list(range(1, len(sequence)+1))
        fig.add_trace(go.Scatter(
            x=positions, 
            y=[0]*len(sequence),
            mode='markers+text',
            text=list(sequence),
            textposition="top center",
            marker=dict(size=15, color='lightgrey')
        ))
        
        # Highlight primer binding sites
        for i, primer in enumerate(primers):
            start = primer['start']
            end = start + primer['length']
            fig.add_trace(go.Scatter(
                x=list(range(start, end)),
                y=[0.5]*(end-start),
                mode='lines+text',
                text=[primer['sequence']],
                line=dict(width=10, color='blue' if i==0 else 'red'),
                name=f"{'Forward' if i==0 else 'Reverse'} Primer"
            ))
        
        fig.update_layout(
            title="Primer Binding Sites",
            showlegend=True,
            yaxis=dict(visible=False),
            height=300,
            hovermode="x unified"
        )
        return fig.to_json()

    def plot_melting_profile(self, sequence, window_size=50):
        """Plot Tm across sequence windows"""
        windows = [sequence[i:i+window_size] for i in range(0, len(sequence)-window_size)]
        tm_values = [self.analyzer.calculate_melting_temperature(w) for w in windows]
        gc_content = [calculate_gc_content(w) for w in windows]
        
        fig = make_subplots(specs=[[{"secondary_y": True}]])
        fig.add_trace(go.Scatter(
            y=tm_values, 
            name="Melting Temp (°C)",
            line=dict(color='red')
        ))
        fig.add_trace(go.Scatter(
            y=gc_content,
            name="GC Content (%)",
            line=dict(color='blue')
        ), secondary_y=True)
        
        fig.update_layout(
            title="Melting Temperature Profile",
            xaxis_title="Sequence Position",
            yaxis_title="Temperature (°C)",
            height=400
        )
        return fig.to_json()

    def plot_amplicon_distribution(self, amplicons):
        """Histogram of amplicon sizes"""
        sizes = [a['amplicon_size'] for a in amplicons]
        fig = go.Figure(data=[go.Histogram(x=sizes, nbinsx=20)])
        fig.update_layout(
            title="Amplicon Size Distribution",
            xaxis_title="Product Size (bp)",
            yaxis_title="Frequency",
            height=300
        )
        return fig.to_json()

    def visualize_dimer(self, dimer_data):
        """Generate dimer structure visualization"""
        fig = go.Figure()
        
        # Add base pairing visualization
        for i, (base1, base2) in enumerate(zip(dimer_data['primer1_region'], 
                                              dimer_data['primer2_region'])):
            fig.add_trace(go.Scatter(
                x=[i, i],
                y=[0, 1],
                mode='lines+markers+text',
                text=[base1, base2],
                textposition="top center",
                line=dict(width=2, color='green' if base1==base2 else 'red'),
                marker=dict(size=15)
            ))
        
        fig.update_layout(
            title=f"Primer Dimer Structure (Score: {dimer_data['stability_score']})",
            showlegend=False,
            yaxis=dict(visible=False),
            height=200
        )
=======
import plotly.graph_objects as go
from plotly.subplots import make_subplots

def calculate_gc_content(sequence):
    """Calculate GC content percentage of a sequence"""
    gc_count = sequence.count('G') + sequence.count('C')
    return (gc_count / len(sequence)) * 100 if sequence else 0

from .algorithms import PrimerAnalyzer

class VisualizationEngine:
    def __init__(self):
        self.analyzer = PrimerAnalyzer()

    def plot_binding_sites(self, sequence, primers):
        """Visualize primer binding positions on target sequence"""
        fig = go.Figure()
        
        # Add sequence as x-axis
        positions = list(range(1, len(sequence)+1))
        fig.add_trace(go.Scatter(
            x=positions, 
            y=[0]*len(sequence),
            mode='markers+text',
            text=list(sequence),
            textposition="top center",
            marker=dict(size=15, color='lightgrey')
        ))
        
        # Highlight primer binding sites
        for i, primer in enumerate(primers):
            start = primer['start']
            end = start + primer['length']
            fig.add_trace(go.Scatter(
                x=list(range(start, end)),
                y=[0.5]*(end-start),
                mode='lines+text',
                text=[primer['sequence']],
                line=dict(width=10, color='blue' if i==0 else 'red'),
                name=f"{'Forward' if i==0 else 'Reverse'} Primer"
            ))
        
        fig.update_layout(
            title="Primer Binding Sites",
            showlegend=True,
            yaxis=dict(visible=False),
            height=300,
            hovermode="x unified"
        )
        return fig.to_json()

    def plot_melting_profile(self, sequence, window_size=50):
        """Plot Tm across sequence windows"""
        windows = [sequence[i:i+window_size] for i in range(0, len(sequence)-window_size)]
        tm_values = [self.analyzer.calculate_melting_temperature(w) for w in windows]
        gc_content = [calculate_gc_content(w) for w in windows]
        
        fig = make_subplots(specs=[[{"secondary_y": True}]])
        fig.add_trace(go.Scatter(
            y=tm_values, 
            name="Melting Temp (°C)",
            line=dict(color='red')
        ))
        fig.add_trace(go.Scatter(
            y=gc_content,
            name="GC Content (%)",
            line=dict(color='blue')
        ), secondary_y=True)
        
        fig.update_layout(
            title="Melting Temperature Profile",
            xaxis_title="Sequence Position",
            yaxis_title="Temperature (°C)",
            height=400
        )
        return fig.to_json()

    def plot_amplicon_distribution(self, amplicons):
        """Histogram of amplicon sizes"""
        sizes = [a['amplicon_size'] for a in amplicons]
        fig = go.Figure(data=[go.Histogram(x=sizes, nbinsx=20)])
        fig.update_layout(
            title="Amplicon Size Distribution",
            xaxis_title="Product Size (bp)",
            yaxis_title="Frequency",
            height=300
        )
        return fig.to_json()

    def visualize_dimer(self, dimer_data):
        """Generate dimer structure visualization"""
        fig = go.Figure()
        
        # Add base pairing visualization
        for i, (base1, base2) in enumerate(zip(dimer_data['primer1_region'], 
                                              dimer_data['primer2_region'])):
            fig.add_trace(go.Scatter(
                x=[i, i],
                y=[0, 1],
                mode='lines+markers+text',
                text=[base1, base2],
                textposition="top center",
                line=dict(width=2, color='green' if base1==base2 else 'red'),
                marker=dict(size=15)
            ))
        
        fig.update_layout(
            title=f"Primer Dimer Structure (Score: {dimer_data['stability_score']})",
            showlegend=False,
            yaxis=dict(visible=False),
            height=200
        )
>>>>>>> 75f1c6a80d6b193bf8a1fa04c4fb03060d16c47f
        return fig.to_json()