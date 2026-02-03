# Computer Vision & Pose Estimation Skill

**RTMPose + OpenVINO** für digitalTwin - Ergonomische Arbeitsplatzanalyse

## Wann aktiviert?

- Keywords: RTMPose, OpenVINO, pose estimation, REBA, RULA, ergonomic, keypoint, OpenCV
- Du arbeitest mit `models/`, Computer Vision Code
- Pose Estimation, Ergonomie-Analyse

## System-Überblick

```
Camera Stream → RTMPose (ONNX) → 17 Keypoints → REBA/RULA Score
                    ↓
               OpenVINO Runtime (CPU/GPU Optimization)
```

---

## 1. RTMPose Model Setup

### Model Download

```python
# models/download_model.py
import urllib.request
import os

MODEL_URL = "https://download.openmmlab.com/mmpose/v1/projects/rtmposev1/onnx_sdk/rtmpose-m_simcc-body7_pt-body7_420e-256x192-e48f03d0_20230504.onnx"
MODEL_PATH = "models/rtmpose-m.onnx"

def download_model():
    if not os.path.exists(MODEL_PATH):
        print(f"Downloading RTMPose model...")
        urllib.request.urlretrieve(MODEL_URL, MODEL_PATH)
        print(f"Model saved to {MODEL_PATH}")
    else:
        print("Model already exists")

if __name__ == "__main__":
    download_model()
```

```bash
python models/download_model.py
```

### Model Specs

- **Eingabe:** 256×192 RGB Image
- **Ausgabe:** 17 Keypoints (COCO Format)
- **Performance:** 90+ FPS (Intel i7-11700 CPU)

**17 Keypoints (COCO):**
0. Nose
1. Left Eye
2. Right Eye
3. Left Ear
4. Right Ear
5. Left Shoulder
6. Right Shoulder
7. Left Elbow
8. Right Elbow
9. Left Wrist
10. Right Wrist
11. Left Hip
12. Right Hip
13. Left Knee
14. Right Knee
15. Left Ankle
16. Right Ankle

---

## 2. OpenVINO Integration

### Setup

```bash
pip install openvino openvino-dev opencv-python numpy
```

### Model Loading

```python
# backend/services/pose_service.py
from openvino.runtime import Core
import numpy as np
import cv2

class PoseEstimator:
    def __init__(self, model_path: str):
        self.core = Core()

        # Load model
        self.model = self.core.read_model(model_path)
        self.compiled_model = self.core.compile_model(self.model, "CPU")

        # Get I/O info
        self.input_layer = self.compiled_model.input(0)
        self.output_layer = self.compiled_model.output(0)

        # Input shape: [1, 3, 256, 192]
        self.input_height = 256
        self.input_width = 192

    def preprocess(self, frame: np.ndarray) -> np.ndarray:
        """Preprocess image for RTMPose"""
        # Resize to model input size
        resized = cv2.resize(frame, (self.input_width, self.input_height))

        # Normalize to [0, 1]
        normalized = resized.astype(np.float32) / 255.0

        # HWC → CHW (Height, Width, Channels → Channels, Height, Width)
        transposed = np.transpose(normalized, (2, 0, 1))

        # Add batch dimension: [3, 256, 192] → [1, 3, 256, 192]
        batched = np.expand_dims(transposed, axis=0)

        return batched

    def infer(self, frame: np.ndarray) -> np.ndarray:
        """Run inference"""
        input_data = self.preprocess(frame)

        # Inference
        result = self.compiled_model([input_data])[self.output_layer]

        return result

    def postprocess(self, output: np.ndarray, original_shape: tuple) -> list:
        """Convert model output to keypoints"""
        keypoints = []

        # RTMPose Output Shape: [1, 17, 64, 48] (heatmaps)
        # Extract max positions from heatmaps
        for i in range(17):
            heatmap = output[0, i, :, :]

            # Find max position
            max_idx = np.unravel_index(np.argmax(heatmap), heatmap.shape)
            y, x = max_idx

            # Scale back to original image size
            orig_h, orig_w = original_shape
            x_scaled = x * (orig_w / 48)  # Heatmap width
            y_scaled = y * (orig_h / 64)  # Heatmap height

            confidence = heatmap[y, x]

            keypoints.append({
                'id': i,
                'x': float(x_scaled),
                'y': float(y_scaled),
                'confidence': float(confidence)
            })

        return keypoints
```

---

## 3. REBA/RULA Scoring

### REBA (Rapid Entire Body Assessment)

**Bewertet:** Ganzkörper-Haltung (stehend/sitzend)

```python
def calculate_reba_score(keypoints: list) -> dict:
    """
    REBA Score: 1 (niedrig) - 15 (sehr hoch)
    """
    # Extract key positions
    left_shoulder = keypoints[5]
    right_shoulder = keypoints[6]
    left_elbow = keypoints[7]
    right_elbow = keypoints[8]
    left_wrist = keypoints[9]
    right_wrist = keypoints[10]
    left_hip = keypoints[11]
    right_hip = keypoints[12]

    # Calculate angles
    neck_angle = calculate_angle(left_shoulder, right_shoulder, [left_shoulder['x'], 0])
    trunk_angle = calculate_trunk_angle(left_shoulder, left_hip)
    leg_score = assess_leg_position(left_hip, keypoints[13])  # Left knee

    # REBA Table A (Neck, Trunk, Legs)
    table_a_score = get_table_a_score(neck_angle, trunk_angle, leg_score)

    # REBA Table B (Arms, Wrists)
    upper_arm_angle = calculate_angle(left_shoulder, left_elbow, left_hip)
    lower_arm_angle = calculate_angle(left_elbow, left_wrist, left_shoulder)
    wrist_angle = calculate_wrist_angle(left_elbow, left_wrist)

    table_b_score = get_table_b_score(upper_arm_angle, lower_arm_angle, wrist_angle)

    # Combine scores
    reba_score = combine_reba_scores(table_a_score, table_b_score)

    return {
        'score': reba_score,
        'risk_level': get_risk_level(reba_score),
        'details': {
            'neck': neck_angle,
            'trunk': trunk_angle,
            'upper_arm': upper_arm_angle
        }
    }

def get_risk_level(score: int) -> str:
    if score <= 1:
        return 'Negligible'
    elif score <= 3:
        return 'Low'
    elif score <= 7:
        return 'Medium'
    elif score <= 10:
        return 'High'
    else:
        return 'Very High'

def calculate_angle(p1: dict, p2: dict, p3: dict) -> float:
    """Calculate angle between 3 points"""
    v1 = np.array([p1['x'] - p2['x'], p1['y'] - p2['y']])
    v2 = np.array([p3['x'] - p2['x'], p3['y'] - p2['y']])

    cos_angle = np.dot(v1, v2) / (np.linalg.norm(v1) * np.linalg.norm(v2))
    angle = np.degrees(np.arccos(np.clip(cos_angle, -1.0, 1.0)))

    return angle
```

---

## 4. FastAPI Integration

```python
# backend/main.py
from fastapi import FastAPI, File, UploadFile
from fastapi.responses import JSONResponse
import cv2
import numpy as np
from services.pose_service import PoseEstimator, calculate_reba_score

app = FastAPI()

# Initialize model
pose_estimator = PoseEstimator("models/rtmpose-m.onnx")

@app.post("/api/analyze-pose")
async function analyze_pose(file: UploadFile = File(...)):
    # Read image
    contents = await file.read()
    nparr = np.frombuffer(contents, np.uint8)
    frame = cv2.imdecode(nparr, cv2.IMREAD_COLOR)

    # Estimate pose
    output = pose_estimator.infer(frame)
    keypoints = pose_estimator.postprocess(output, frame.shape[:2])

    # Calculate REBA
    reba = calculate_reba_score(keypoints)

    return JSONResponse({
        'keypoints': keypoints,
        'reba_score': reba['score'],
        'risk_level': reba['risk_level'],
        'details': reba['details']
    })
```

---

## 5. Frontend Integration (React)

```tsx
// components/CameraView.tsx
import { useRef, useState } from 'react';

export function CameraView() {
  const videoRef = useRef<HTMLVideoElement>(null);
  const canvasRef = useRef<HTMLCanvasElement>(null);
  const [result, setResult] = useState<any>(null);

  async function captureAndAnalyze() {
    const video = videoRef.current;
    const canvas = canvasRef.current;
    if (!video || !canvas) return;

    // Capture frame
    const ctx = canvas.getContext('2d')!;
    canvas.width = video.videoWidth;
    canvas.height = video.videoHeight;
    ctx.drawImage(video, 0, 0);

    // Convert to Blob
    const blob = await new Promise<Blob>((resolve) =>
      canvas.toBlob(resolve as BlobCallback, 'image/jpeg')
    );

    // Send to API
    const formData = new FormData();
    formData.append('file', blob!);

    const response = await fetch('/api/analyze-pose', {
      method: 'POST',
      body: formData
    });

    const data = await response.json();
    setResult(data);
  }

  return (
    <div>
      <video ref={videoRef} autoPlay />
      <canvas ref={canvasRef} style={{ display: 'none' }} />

      <button onClick={captureAndAnalyze}>Analyze Posture</button>

      {result && (
        <div>
          <h3>REBA Score: {result.reba_score}</h3>
          <p>Risk Level: {result.risk_level}</p>
        </div>
      )}
    </div>
  );
}
```

---

## 6. Performance Optimization

### GPU Acceleration

```python
# Use GPU wenn verfügbar
self.compiled_model = self.core.compile_model(self.model, "GPU")  # oder "CPU"
```

### Batch Processing

```python
def infer_batch(self, frames: list[np.ndarray]) -> list[np.ndarray]:
    """Process multiple frames at once"""
    batch = np.stack([self.preprocess(f) for f in frames])
    result = self.compiled_model([batch])[self.output_layer]
    return result
```

### Frame Skipping

```tsx
// Analysiere nur jedes N-te Frame
let frameCount = 0;

function onFrame() {
  frameCount++;
  if (frameCount % 30 === 0) {  // Alle 30 Frames (1x/Sekunde bei 30fps)
    captureAndAnalyze();
  }
}
```

---

## 7. Mock-Daten für Development

```python
# backend/services/mock_pose_service.py
def get_mock_keypoints() -> list:
    """Mock keypoints für Development ohne Camera"""
    return [
        {'id': 0, 'x': 128, 'y': 50, 'confidence': 0.95},   # Nose
        {'id': 5, 'x': 100, 'y': 80, 'confidence': 0.92},   # Left Shoulder
        {'id': 6, 'x': 156, 'y': 80, 'confidence': 0.91},   # Right Shoulder
        # ... 14 weitere
    ]
```

---

## Troubleshooting

**Model lädt nicht?**
→ Prüfe Pfad, lade Model erneut herunter

**Niedrige FPS?**
→ Nutze OpenVINO GPU, reduziere Input-Auflösung

**Falsche Keypoints?**
→ Prüfe Preprocessing (Normalization, Transpose)

## Ressourcen

- [RTMPose Paper](https://arxiv.org/abs/2303.07399)
- [OpenVINO Docs](https://docs.openvino.ai/)
- [REBA Guide](https://ergo-plus.com/reba-assessment-tool-guide/)
